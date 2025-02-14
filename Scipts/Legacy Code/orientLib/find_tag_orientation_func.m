function [rot_mat, in_pct] = find_tag_orientation_func(Acc, Depth, sample_freq,...
  surf_depth, plot_interval, method)
% Find the rotation matrix 'rot_mat' that rotates the xyz axes of tag  
% sensor readings (accelerometer, magmetometer, gyroscope), so that the 
% rotated sensor axes are aligned with the animal's body frame: x-surge  
% points forward, y-sway points left, z-heave points up (of the animal). 
% Be aware that in this coordinate system, the positive rotation around 
% y-sway axis is acutally pitching down.
%
% The method works under the assumptions that:
% 1) The tag is NOT sliding on the animal for the provided section of data.
%    i.e. the relative tag oriantation is unkonwn but FIXED.
% 2) For the *majority* of the section of data, the animal is swimming
%    normally, rather than rolling to the left or right.
% 3) The animal is cetacea (e.g. dolphin, whale etc.) or demonstrates 
%    very similar movement/gait pattern.
%
% The method works for the case the animal rolling to the left or right, as
% long as those are not dominating the entire data section. That means, 
% the method works better when: 1) the data section is long and 2) the 
% animal is swimming 'normally' rather than a strange palyful way.
%
% INPUT:
% Acc           n-by-3 tag accelerometer data. [unit does not matter]
% Depth         n-by-1 *calibrated* tag depth data. [m] 
% sample_freq   1-by-1 scalar, sampleing frequency of 'Acc' and 'Depth'.
% surf_depth    [Optional] 1-by-1 scalar, defines the depth threshold for  
%               animal surfacing events, and excludes those for orientation 
%               anlysis. This is because animal can sometimes perform  
%               confusing non-swimming behavior at the surface, and 
%               confuses the method. 
%               Set to '[]' (i.e. surf_depth = []) to include all data. 
%               The default value is '[]'.
% plot_interval [Optional] 1-by-1 scalar, when visualizing data, data 
%               points are plotted every 'plot_interval', so the plot is 
%               not too dense. I.e. plot(Acc(1:plot_interval:n,:))
%               Set to '[]' (i.e. plot_interval = [], also is default) to 
%               use automatic plot interval, so a fix number of points are 
%               shown. Set to 0 to turn off plot.
%               *Set to '0' (i.e. plot_interval = 0) to turn off plots.*
% method        [Optional] 1-by-1 scalar value: '1', '2' or '[]' indicates 
%               the method to use when fitting to data. "1" works for the
%               case that the animal rolls more and doing less deep dive.
%               '2' works for the case that the animal does more deep dive.
%               '[]' lets the program to decide which one to use by
%               checking which one gives more fitting inliers.
%               The default value is '[]', letting the program run in auto.
%               Check the support function "fit_fad" in this script for 
%               more detail. 
% OUTPUT:
% rot_mat       3-by-3 rotation matrix, that rotates tag-frame sensor data
%               into animal-frame as described above. Example usage:
%               Acc_rotated = Acc*rot_mat;
%               Mag_rotated = Mag*rot_mat;
% in_pct        1-by-1 scalar, model fitting inlier percentage, ranges [0,1]. 
%               Higher percentage indicates better fits. Note that due to 
%               how RANSAC fitting is used here, certain amount of outliers
%               are normal, a good fit does not necessarily means a very 
%               high "in_pct" value.
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 07/07/2020
% 04/05/2022: Added a line to low pass filter depth.
% ====================


%% Manage input.
% Check if exist nan in Acc or Depth.
assert(sum(sum(isnan(Acc))) == 0 & sum(isnan(Depth)) == 0,...
  "Exist NaN in the input accelerometer and/or depth data, please fix input.")

% Default values.
if nargin == 3
  surf_depth = [];
  plot_interval = [];
  method = [];
elseif nargin == 4
  plot_interval = [];
  method = [];
elseif nargin == 5
  method = [];
end


% Check if method is valid.
assert(isempty(method) || method == 1 || method == 2,...
  'Invalid method specification.')

% Automatic plot interval.
if isempty(plot_interval)
  n_points_show = 10000;
  plot_interval = ceil(length(Acc)/n_points_show);
end

% Plot on/off switch.
plot_on = true;
if plot_interval <= 0
  plot_on = false;
else
  plot_interval = ceil(plot_interval);
end


%% Prepare data.
% Flip Depth if they are positive. In this script, all depth (i.e. z-axis
% position) are assumed negative or zero, so when the animal descends, depth
% decreases.
if nanmean(Depth) > 0
  Depth = -Depth;
end

% Filter depth.
Depth = movmean(Depth, ceil(0.5*sample_freq));

% Depth speed.
Depth_speed = diff(Depth)*sample_freq;
Depth_speed(end+1) = Depth_speed(end);

% Flip input surface threshold if needed, to align with depth direction.
% Surface threshold is used to exclude surface data.
if surf_depth > 0
  surf_depth = -surf_depth;
end
if isempty(surf_depth)
  surf_depth = Inf; % Surface inf high <==> Include all data.
end

% Filter and scale down Acc to unit [g].
A0 = movmean(Acc, ceil(0.5*sample_freq));
A_norm = sqrt(sum(A0.^2, 2));
A0 = A0./mean(A_norm);

% Animal flat/descend/ascend definitions.
flat_thrs = max(0.05, prctile(Depth_speed, 70));
% flat_thrs = max(0.02, prctile(Depth_speed, 30));

ascend_thrs = max(0.08, prctile(Depth_speed, 80));
% ascend_thrs = max(0.05, prctile(Depth_speed, 70));

idx_flat = abs(Depth_speed) <= flat_thrs; % prctile(Depth_speed, 70)
idx_descend = Depth_speed < -ascend_thrs; % prctile(Depth_speed, 12)
idx_ascend = Depth_speed > ascend_thrs; % prctile(Depth_speed, 88)

% Other concerns, excluding surfaces.
idx_other = Depth < surf_depth;


%% Find the dominant Acc (i.e. gravity) direction
% corresponding to the animal swimming flat/ascend/descend.
% By fitting two planes to Acc data corresponding to animal swimming 
% flat/ascend/descend. Using RANSAC fitting techniques, which is robust to
% outliers.
if isempty(method)
  % Try both fitting methods, use results from high "in_pct" one.
  m = 1;
  [Af0_temp{m}, Aa0_temp{m}, Ad0_temp{m}, model_1_temp{m},model_2_temp{m},...
    in_pct_temp{m}] = fit_fad(m);
  m = 2;
  [Af0_temp{m}, Aa0_temp{m}, Ad0_temp{m}, model_1_temp{m}, model_2_temp{m},...
    in_pct_temp{m}] = fit_fad(m);
  if in_pct_temp{1} > in_pct_temp{2}
    method = 1;
  else
    method = 2;
  end
  Af0 = Af0_temp{method};
  Aa0 = Aa0_temp{method};
  Ad0 = Ad0_temp{method};
  model_1 = model_1_temp{method};
  model_2 = model_2_temp{method};
  in_pct = in_pct_temp{method};
else
  [Af0, Aa0, Ad0, model_1, model_2, in_pct] = fit_fad(method);
end

  
%% Perform rotations to align data. A0 -> A1 (z-axis align) -> A2 (all align).
% Find (normalized) dominant gravity direction when the animal is swimming 
% flat/ascend/descend.
scale_val = 1.35;
A0_flat_g = scale_val * Af0./norm(Af0);    %mean(A0(idx_flat,:)); % Naive way of doing this.
A0_ascend_g = scale_val * Aa0./norm(Aa0);  %mean(A0(idx_ascend,:));
A0_descend_g = scale_val * Ad0./norm(Ad0); %mean(A0(idx_descend,:));

% First rotations to find A1. 
% Find rotation to align A_g with the vertical axis of tag (i.e. z-axis).
% and rotate data accordingly.
A_g_tag = [0, 0, 1];
rot_mat_1 = vrrotvec2mat_local(vrrotvec_local(A_g_tag, A0_flat_g));

% Second rotation to find A2.
% Now vertical (z) aixs direction is fixed, let's fix the rotation around it.
% Use the section of data corresponding to animal descend.
% Find rotation around z so when descending, gravity direction is aligned
% with negative direction of x-axis.
A1_descend_g_xy = A0_descend_g*rot_mat_1;
A1_descend_g_xy(1,3) = 0;

A1_descend_g_tag = [-1, 0, 0];

% Overall rotation matrix.

% If A1_descend_g_xy has NaNs, discard and use previous rotation matrix
if any(isnan(A1_descend_g_xy))
    rot_mat = [];
    % Show failure to compute rotation.
    disp('Tag orientation not computed, using previous rotation.')
else
    rot_mat_2 = vrrotvec2mat_local(vrrotvec_local(A1_descend_g_tag, A1_descend_g_xy));
    rot_mat = rot_mat_1*rot_mat_2;
    % Show method used and success.
    disp(['Tag orientation found, method used: ', num2str(method)])
end

% Use format. With X_raw (n-by-3) be the sensor axis needs to be fixed.
% X_fix = X_raw*rot_mat;

%% Visualizations. 
if plot_on
    if isempty(rot_mat)
       disp('Secondary plots not run: section orientation not found.')
    else
        plot_all()
    end
end

%% Sub-functions. =========================================================
function [Af0, Aa0, Ad0, model_1, model_2, in_pct] = fit_fad(method)
  % Two methods to find the dominant Acc (i.e. gravity) directions
  % corresponding to animal swimming flat ("Af0"), ascend ("Aa0") and 
  % descend ("Ad0"). By fitting two planes to data ("model_1" and
  % "model_2"). Fitting inlier percentage is return as "in_pct", higher
  % inlier percentage means better fit.
  %
  % Method 1 - fit flat first: Fit a plane to flat (f) data first, then fit
  % a second plane to ascend and descend (ad) data, the second plane is also 
  % perpendicular to the first plane.
  %
  % Method 2 - fit ad first: Fit a plane to ascend and descend (ad) data 
  % first, then fit a second plane to flat (f) data, the second plane is 
  % also perpendicular to the first plane.
  if nargin == 0
    method = 1; % Fit to flat first as default.
  end
  assert(method == 1 | method == 2, 'Invalid method specification.')
  
  if method == 1
    idx_analysis_1 = idx_flat&idx_other;
    idx_analysis_2 = (idx_flat|idx_ascend|idx_descend)&idx_other;  
  else
    %idx_analysis_1 = (idx_ascend|idx_descend)&idx_other;
    idx_analysis_1 = (idx_flat|idx_ascend|idx_descend)&idx_other;
    idx_analysis_2 = idx_flat&idx_other;
  end

  % RANSAC plane fitting, the planes go through zero (0,0,0).
  % Fit first plane.
  data_1 = A0(idx_analysis_1, :);
  fitFcn_1 = @(points) [points(:,1:2)\points(:,3); 0];
  distFcn = @(model, points) abs(points*[model(1:2);-1])./...
    norm([model(1:2); -1]);
  sampleSize_1 = 2;
  maxDistance = 0.1;

  [model_1, idx_in_1] = ransac_local(data_1,fitFcn_1,distFcn,sampleSize_1,...
    maxDistance, 'MaxSamplingAttempts', 1e10, 'MaxNumTrials', 1e10,...
    'Confidence', 99.999999);

  % Fit a second plane.
  norm_1 = [model_1(1:2); -1]';
  norm_1 = norm_1./norm(norm_1);
  fitFcn_2 = @(point) fitPlaneDir(norm_1, point);
  sampleSize_2 = 1;
  data_2 = A0(idx_analysis_2, :);

  [model_2, idx_in_2] = ransac_local(data_2,fitFcn_2,distFcn,sampleSize_2,...
    maxDistance, 'MaxSamplingAttempts', 1e10, 'MaxNumTrials', 1e10,...
    'Confidence', 99.999999);

  % Manage the indices for finding inlier data points.
  idx_num = (1:length(A0))';
  idx_num_1 = idx_num(idx_analysis_1);
  idx_num_1_in = idx_num_1(idx_in_1);

  idx_num_2 = idx_num(idx_analysis_2);
  idx_num_2_in = idx_num_2(idx_in_2);
  
  idx_num_12 = union(idx_num_1, idx_num_2);
  idx_num_12_in = union(idx_num_1_in, idx_num_2_in);
  
  % Plane fitting inlier percentage.
  in_pct = length(idx_num_12_in)/length(idx_num_12);

  % Find domninant direction using fitted plane.
  A0_flat_mean = mean(A0(idx_flat&idx_other,:));
  A0_descend_mean = mean(A0(idx_descend&idx_other,:));
  A0_ascend_mean = mean(A0(idx_ascend&idx_other,:));
  norm_2 = [model_2(1:2); -1]';
  norm_2 = norm_2./norm(norm_2);
  
  % Project average A0 values of flat/ascend/descend to planes, to find the
  % dominant directions of A0 corresponding to animal swimming flat/ascend
  % /descend.
  Af0 = projectVectorToPlane(norm_1, A0_flat_mean);
  Af0 = projectVectorToPlane(norm_2, Af0);
  
  if method == 1
    Ad0 = projectVectorToPlane(norm_2, A0_descend_mean);
    Aa0 = projectVectorToPlane(norm_2, A0_ascend_mean);    
  else
    Ad0 = projectVectorToPlane(norm_1, A0_descend_mean);
    Aa0 = projectVectorToPlane(norm_1, A0_ascend_mean);
  end
end


function plot_all()
  % Visualizations.
  %% Plot supporting computations.
  norm_1 = [model_1(1:2); -1]';
  norm_1 = norm_1./norm(norm_1);
  norm_2 = [model_2(1:2); -1]';
  norm_2 = norm_2./norm(norm_2);
  
  % Mesh of plane.
  acc_bound = 1.2;
  [X_mesh,Y_mesh] = meshgrid(linspace(-acc_bound, acc_bound,50),...
    linspace(-acc_bound, acc_bound, 50));
  Z_mesh_1 = model_1(1)*X_mesh + model_1(2)*Y_mesh;
  Z_mesh_1(abs(Z_mesh_1) > acc_bound) = nan;

  Z_mesh_2 = model_2(1)*X_mesh + model_2(2)*Y_mesh;
  Z_mesh_2(abs(Z_mesh_2) > acc_bound) = nan;
  
  % Dominant directions.
  A1_flat_g = A0_flat_g*rot_mat_1;
  A1_descend_g = A0_descend_g*rot_mat_1;
  A1_ascend_g = A0_ascend_g*rot_mat_1;
  
  A2_flat_g = A1_flat_g*rot_mat_2;
  A2_descend_g = A1_descend_g*rot_mat_2;
  A2_ascend_g = A1_ascend_g*rot_mat_2;
  
  % A1, A2 data.
  A1 = A0*rot_mat_1;
  A2 = A1*rot_mat_2;
  
  % Ordinary least square fit to data.
  if method == 1
    data_ls = A0(idx_flat&idx_other,:);
  else
    data_ls = A0((idx_ascend|idx_descend)&idx_other,:);
  end
    
  model_ols = data_ls(:,1:2)\data_ls(:,3); % Estimate Parameters
  Z_mesh_ols = model_ols(1)*X_mesh + model_ols(2)*Y_mesh;
  Z_mesh_ols(abs(Z_mesh_ols) > acc_bound) = nan;

  
  %% Plot A0 raw.
  idx_plot = false(length(A0), 1);
  idx_plot(1:plot_interval:length(A0)) = true;

  figure
  plot3(A0(idx_plot, 1),A0(idx_plot, 2),...
    A0(idx_plot, 3),'.')

  grid on
  axis equal
  xlabel('x-surge')
  ylabel('y-sway')
  zlabel('z-heave')
  title('A0, raw')
  legend('A0-all')

  %% Plot A0 colored.
  figure
  plot3(A0(idx_plot & idx_flat & idx_other, 1),...
        A0(idx_plot & idx_flat & idx_other, 2),...
        A0(idx_plot & idx_flat & idx_other, 3),'.')
  hold on
  plot3(A0(idx_plot & idx_descend & idx_other,1),...
        A0(idx_plot & idx_descend & idx_other,2),...
        A0(idx_plot & idx_descend & idx_other,3),'.')
  plot3(A0(idx_plot & idx_ascend & idx_other,1),...
        A0(idx_plot & idx_ascend & idx_other,2),...
        A0(idx_plot & idx_ascend & idx_other,3),'.')
  plot3([0, A0_flat_g(1,1)], [0, A0_flat_g(1,2)], [0, A0_flat_g(1,3)],...
    '-*', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A0_descend_g(1,1)], [0, A0_descend_g(1,2)], [0, A0_descend_g(1,3)],...
    '-.o', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A0_ascend_g(1,1)], [0, A0_ascend_g(1,2)], [0, A0_ascend_g(1,3)],...
    '-->', 'lineWidth',2, 'markerSize', 8)

  grid on
  axis equal
  xlabel('x-surge')
  ylabel('y-sway')
  zlabel('z-heave')
  title('A0, colored')
  legend('A0-flat', 'A0-descend', 'A0-ascend', 'A0-flat-g', 'A0-descend-g',...
      'A0-ascend-g')
    
  %% Plot Least square vs. Ransac fit to flat data.
  figure
  plot3(data_ls(1:plot_interval:length(data_ls),1),...
        data_ls(1:plot_interval:length(data_ls),2),...
        data_ls(1:plot_interval:length(data_ls),3),'.')
  hold on

  pl_ls = meshc(X_mesh, Y_mesh, Z_mesh_ols);
  pl_ls(1).EdgeColor = 'c';
  pl_ls(2).EdgeColor = 'w';

  pl_ransac = meshc(X_mesh, Y_mesh, Z_mesh_1);
  pl_ransac(1).EdgeColor = 'm';
  pl_ransac(2).EdgeColor = 'w';

  hold off
  grid on
  xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)');
  title('Ordinary Least Square vs. RANSAC');
  grid on
  axis equal
  legend('A0', 'Ordinary least square fit', '', 'RANSAC fit', '')

  %% Plot A0 plane fit.
  figure
  plot3(A0(idx_plot & idx_flat & idx_other,1),...
        A0(idx_plot & idx_flat & idx_other,2),...
        A0(idx_plot & idx_flat & idx_other,3),'.')
  hold on
  plot3(A0(idx_plot & idx_descend & idx_other,1),...
        A0(idx_plot & idx_descend & idx_other,2),...
        A0(idx_plot & idx_descend & idx_other,3),'.')
  plot3(A0(idx_plot & idx_ascend & idx_other,1),...
        A0(idx_plot & idx_ascend & idx_other,2),...
        A0(idx_plot & idx_ascend & idx_other,3),'.')
  plot3([0, A0_flat_g(1,1)], [0, A0_flat_g(1,2)], [0, A0_flat_g(1,3)],...
    '-*', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A0_descend_g(1,1)], [0, A0_descend_g(1,2)], [0, A0_descend_g(1,3)],...
    '-.o', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A0_ascend_g(1,1)], [0, A0_ascend_g(1,2)], [0, A0_ascend_g(1,3)],...
    '-->', 'lineWidth',2, 'markerSize', 8)

  pl_ransac = meshc(X_mesh, Y_mesh, Z_mesh_1);
  pl_ransac(1).EdgeColor = 'm';
  pl_ransac(2).EdgeColor = 'w';

  pl_ransac_ad = meshc(X_mesh, Y_mesh, Z_mesh_2);
  pl_ransac_ad(1).EdgeColor = 'r';
  pl_ransac_ad(2).EdgeColor = 'w';

  
  plot3([0, norm_1(1,1)], [0, norm_1(1,2)], [0, norm_1(1,3)],...
    'm-^', 'lineWidth',2, 'markerSize', 8)
  plot3([0, norm_2(1,1)], [0, norm_2(1,2)], [0, norm_2(1,3)],...
    'r-^', 'lineWidth',2, 'markerSize', 8)
  
  grid on
  axis equal
  xlabel('x-surge')
  ylabel('y-sway')
  zlabel('z-heave')
  title('A0, plane fit')
  legend('A0-flat', 'A0-descend', 'A0-ascend', 'A0-flat-g', 'A0-descend-g',...
      'A0-ascend-g', 'Plane-1','',...
      'Plane-2','', 'Normal-1', 'Normal-2')

  %% Plot A1.
  figure
  plot3(A1(idx_plot & idx_flat & idx_other,1),...
        A1(idx_plot & idx_flat & idx_other,2),...
        A1(idx_plot & idx_flat & idx_other,3),'.')
  hold on
  plot3(A1(idx_plot & idx_descend & idx_other,1),...
        A1(idx_plot & idx_descend & idx_other,2),...
        A1(idx_plot & idx_descend & idx_other,3),'.')
  plot3(A1(idx_plot & idx_ascend & idx_other,1),...
        A1(idx_plot & idx_ascend & idx_other,2),...
        A1(idx_plot & idx_ascend & idx_other,3),'.')

  plot3([0, A1_flat_g(1,1)], [0, A1_flat_g(1,2)], [0, A1_flat_g(1,3)],...
    '-*', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A1_descend_g(1,1)], [0, A1_descend_g(1,2)], [0, A1_descend_g(1,3)],...
    '-.o', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A1_ascend_g(1,1)], [0, A1_ascend_g(1,2)], [0, A1_ascend_g(1,3)],...
    '-->', 'lineWidth',2, 'markerSize', 8)

  grid on
  axis equal
  xlabel('x-surge')
  ylabel('y-sway')
  zlabel('z-heave')
  title('A1, z-axis aligned')
  legend('A1-flat', 'A1-descend', 'A1-ascend', 'A1-flat-g','A1-descend-g',...
      'A1-ascend-g')

  %% Plot A2.
  % Number of points plotted. 
  n_flat = sum(idx_plot & idx_flat & idx_other);
  n_asc = sum(idx_plot & idx_ascend & idx_other);
  n_des = sum(idx_plot & idx_descend & idx_other);
  n_fad = n_flat + n_asc + n_des;
  figure
  plot3(A2(idx_plot & idx_flat & idx_other,1),...
        A2(idx_plot & idx_flat & idx_other,2),...
        A2(idx_plot & idx_flat & idx_other,3),'.')
  hold on
  plot3(A2(idx_plot & idx_descend & idx_other,1),...
        A2(idx_plot & idx_descend & idx_other,2),...
        A2(idx_plot & idx_descend & idx_other,3),'.')
  plot3(A2(idx_plot & idx_ascend & idx_other,1),...
        A2(idx_plot & idx_ascend & idx_other,2),...
        A2(idx_plot & idx_ascend & idx_other,3),'.')

  plot3([0, A2_flat_g(1,1)], [0, A2_flat_g(1,2)], [0, A2_flat_g(1,3)],...
    '-*', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A2_descend_g(1,1)], [0, A2_descend_g(1,2)], [0, A2_descend_g(1,3)],...
    '-.o', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A2_ascend_g(1,1)], [0, A2_ascend_g(1,2)], [0, A2_ascend_g(1,3)],...
    '-->', 'lineWidth',2, 'markerSize', 8)
  
  grid on
  axis equal
  xlabel('x-surge')
  ylabel('y-sway')
  zlabel('z-heave')
  title([{'A2, all aligned, number of points showing:'},...
    ['flat (', num2str(n_flat), '), ascend (', num2str(n_asc),...
     '), descend (',num2str(n_des), '), total (', num2str(n_fad),')']])

  legend('A2-flat', 'A2-descend', 'A2-ascend', 'A2-flat-g','A2-descend-g',...
      'A2-ascend-g')

  %% Plot A0 vs. A2.
  idx_plot = false(length(A2), 1);
  idx_plot(1:5*plot_interval:length(A2)) = true;

  figure
  hold on
  plot3(A0(idx_plot,1),A0(idx_plot,2),A0(idx_plot,3),'.')
  plot3(A2(idx_plot,1),A2(idx_plot,2),A2(idx_plot,3),'.')

  plot3([0, A0_flat_g(1,1)], [0, A0_flat_g(1,2)], [0, A0_flat_g(1,3)],...
    'b-*', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A0_descend_g(1,1)], [0, A0_descend_g(1,2)],...
    [0, A0_descend_g(1,3)], 'b-.o', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A0_ascend_g(1,1)], [0, A0_ascend_g(1,2)],...
    [0, A0_ascend_g(1,3)], 'b-->', 'lineWidth',2, 'markerSize', 8)

  plot3([0, A2_flat_g(1,1)], [0, A2_flat_g(1,2)], [0, A2_flat_g(1,3)],...
    'r-*', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A2_descend_g(1,1)], [0, A2_descend_g(1,2)],...
    [0, A2_descend_g(1,3)], 'r-.o', 'lineWidth',2, 'markerSize', 8)
  plot3([0, A2_ascend_g(1,1)], [0, A2_ascend_g(1,2)],...
    [0, A2_ascend_g(1,3)], 'r-->', 'lineWidth',2, 'markerSize', 8)

  grid on
  axis equal
  xlabel('x-surge')
  ylabel('y-sway')
  zlabel('z-heave')
  title('A0 vs. A2')
  legend('A0', 'A2', 'A0-flat-g', 'A0-descend-g', 'A0-ascend-g',...
    'A2-g', 'A2-descend-g', 'A2-ascend-g')
end
end


%% === Helper functions ===========================================
function model = fitPlaneDir(p0, point)
  % Find a plane that must goes through origin and point p0 (and fit with
  % other 'points', or not...).
  % model = [a, b, c], with z = a*x + b*y + c;
  % We know c = 0 already.
  % b = (z0 - a*x0)/y0 
  %   = z0/y0 - a*x0/y0
  % a = (z - z0/y0*y)/(x - x0/y0*y)
  
  normal = cross(p0, point(1,:));
  a = -normal(1)/normal(3);
  b = -normal(2)/normal(3);
  c = 0;
  model = [a; b; c];
end


function v1 = projectVectorToPlane(normal, v0)
  % Project vector 'v0' onto the plane with 'normal', plane goes through 
  % origin. Resulting vector is 'v1'. All 'normal', 'v0', 'v1' start from
  % origin.
  n0 = normal./norm(normal);
  n1 = n0*dot(n0, v0);
  v1 = v0 - n1;
end

