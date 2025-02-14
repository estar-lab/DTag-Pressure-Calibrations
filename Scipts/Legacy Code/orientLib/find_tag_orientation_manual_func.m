function rot_mat = find_tag_orientation_manual_func(Acc, Depth, sample_freq,...
  D_flat, D_ascend, D_descend, surf_depth, plot_interval)
% Same purpose as "find_tag_orientation_manual_func()", but this one asks
% you to define dominant vectors manually. This is for fine scale
% correction of the data.
%
% Find the rotation matrix 'rot_mat' that rotates the xyz axes of tag  
% sensor readings (accelerometer, magmetometer, gyroscope), so that the 
% rotated sensor axes are aligned with the animal's body frame: x-surge  
% points forward, y-sway points left, z-heave points up (of the animal). 
% Be aware that in this coordinate system, the positive rotation around 
% y-sway axis is acutally pitching down.
%
%
% INPUT:
% Acc           n-by-3 tag accelerometer data. [unit does not matter]
% Depth         n-by-1 *calibrated* tag depth data. [m] 
% sample_freq   1-by-1 scalar, sampleing frequency of 'Acc' and 'Depth'.
% D_flat        1-by-3 dominant vector corresponds to flat swim.
% D_ascend      1-by-3 dominant vector corresponds to ascending.
% D_descend     1-by-3 dominant vector corresponds to descending.

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
%               shown. 
%               *Set to '0' (i.e. plot_interval = 0) to turn off plots.*
% OUTPUT:
% rot_mat       3-by-3 rotation matrix, that rotates tag-frame sensor data
%               into animal-frame as described above. Example usage:
%               Acc_rotated = Acc*rot_mat;
%               Mag_rotated = Mag*rot_mat;
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 07/07/2020
% ====================


%% Manage input.
% Check if exist nan in Acc or Depth.
assert(sum(sum(isnan(Acc))) == 0 & sum(isnan(Depth)) == 0,...
  "Exist NaN in the input accelerometer and/or depth data, please fix input.")

% Make sure input dominant vectors have size 1-by-3.
assert(size(D_flat, 1) == 1 && size(D_flat, 2) == 3 &&...
    size(D_ascend, 1) == 1 && size(D_ascend, 2) == 3 &&...
    size(D_descend, 1) == 1 && size(D_descend, 2) == 3,...
    'Input dominant direction vectors need to be 1-by-3 row vectors.')


% Default values.
if nargin == 6
  surf_depth = [];
  plot_interval = [];
elseif nargin == 7
  plot_interval = [];
end

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
flat_thrs = prctile(Depth_speed, 70);
ascend_thrs = prctile(Depth_speed, 80);
idx_flat = abs(Depth_speed) < flat_thrs; % prctile(Depth_speed, 70)
idx_descend = Depth_speed < -ascend_thrs; % prctile(Depth_speed, 12)
idx_ascend = Depth_speed > ascend_thrs; % prctile(Depth_speed, 88)

% Other concerns, excluding surfaces.
idx_other = Depth < surf_depth;

  
%% Perform rotations to align data. A0 -> A1 (z-axis align) -> A2 (all align).
% Find (normalized) dominant gravity direction when the animal is swimming 
% flat/ascend/descend.
scale_val = 1.35;
A0_flat_g = scale_val * D_flat./norm(D_flat);    %mean(A0(idx_flat,:)); % Naive way of doing this.
A0_ascend_g = scale_val * D_ascend./norm(D_ascend);  %mean(A0(idx_ascend,:));
A0_descend_g = scale_val * D_descend./norm(D_descend); %mean(A0(idx_descend,:));

% First rotations to find A1. 
% Find rotation to align A_g with the vertical axis of tag (i.e. z-axis).
% and rotate data accordingly.
A_g_tag = [0, 0, 1];
rot_mat_1 = vrrotvec2mat(vrrotvec(A_g_tag, A0_flat_g));

% Second rotation to find A2.
% Now vertical (z) aixs direction is fixed, let's fix the rotation around it.
% Use the section of data corresponding to animal descend.
% Find rotation around z so when descending, gravity direction is aligned
% with negative direction of x-axis.
A1_descend_g_xy = A0_descend_g*rot_mat_1;
A1_descend_g_xy(1,3) = 0;

A1_descend_g_tag = [-1, 0, 0];
rot_mat_2 = vrrotvec2mat(vrrotvec(A1_descend_g_tag, A1_descend_g_xy));

% Overall rotation matrix.
rot_mat = rot_mat_1*rot_mat_2;
% Use format. With X_raw (n-by-3) be the sensor axis needs to be fixed.
% X_fix = X_raw*rot_mat;


%% Visualizations. 
if plot_on
  plot_all()
end

%% Sub-functions. =========================================================
function plot_all()
  % Visualizations.
  %% Plot supporting computations.
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

end
end
