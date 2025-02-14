function [inliers, fad_status] = find_inliers_func(Acc_1, Depth_1, Acc_2,...
  Depth_2, sample_freq, plot_interval)
% Find the indices of the data points from the 2nd dataset (Acc2, Depth2) 
% that fall within the distribution of the first dataset (Acc1, Depth1).
%
% - Consider cut A1 (and potentially A2) to get denser distribution. 
% Cut flat/ascend/descend separately.
%
%
% INPUT:
% Acc_1         n1-by-3 tag accelerometer data, from the reference 
%               data section. [unit does not matter]
% Depth_1       n1-by-1 calibrated tag depth data, from the reference 
%               data section. [m]
% Acc_2         n2-by-3 tag accelerometer data, from the unknown data
%               section, that to be checked to find the inliers with   
%               respect to the distribution of the reference section. 
%               [unit does not matter]
% Depth_2       n2-by-1 calibrated tag depth data, from the unknown data
%               section. [m]
% sample_freq   1-by-1 scalar, sampleing frequency of 'Acc' and 'Depth'.
% plot_interval 1-by-1 scalar, when visualizing data, data points are
%               plotted every 'plot_interval', so the plot is not too dense.
%               I.e. Acc(1:plot_interval:n,:)
%               Set to '[]' or '0' (i.e. plot_interval = [], or plot_interval 
%               = 0) to turn off plots. The default value is 20.
% OUTPUT:
% inliers       n2-by-1 binary values indicate whether or not a point in
%               "Acc_2" is a inlier in the distribution of "Acc_1".
% fad_status    n2-by-1 scalars indicates whether the animal is flat, 
%               ascend or descend for each Acc_2 point. '1'-flat,
%               '2'-ascend, '3'-descend.
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 07/08/2020
% ====================

%% Manage input.
% Default values.
if nargin == 5
  plot_interval = 20;
end

% Plot on/off switch.
plot_on = true;
if isempty(plot_interval) || plot_interval <= 0
  plot_on = false;
else
  plot_interval = ceil(plot_interval);
end

% Flip Depth if they are positive. In this script, all depth (i.e. z-axis
% position) are negative or zero, so when the animal descends, depth
% decreases.
if nanmean(Depth_1) > 0
  Depth_1 = -Depth_1;
  Depth_2 = -Depth_2;
end


%% Prepare data.
% Filter and scale down Acc to unit [g].
A{1} = movmean(Acc_1, ceil(0.1*sample_freq));
%A{1} = A{1}./mean(sqrt(sum(A{1}.^2, 2)));
A{1} = A{1}./sqrt(sum(A{1}.^2, 2));

A{2} = movmean(Acc_2, ceil(0.1*sample_freq));
%A{2} = A{2}./mean(sqrt(sum(A{2}.^2, 2)));
A{2} = A{2}./sqrt(sum(A{2}.^2, 2));

% Depth speed.
Depth_speed{1} = diff(Depth_1)*sample_freq;
Depth_speed{1}(end+1) = Depth_speed{1}(end);

Depth_speed{2} = diff(Depth_2)*sample_freq;
Depth_speed{2}(end+1) = Depth_speed{2}(end);


% Down sample A{1} and Depth_speed{1} to 10 Hz.
sf_down = 10;
if sample_freq > sf_down
  idx_use = 1:ceil(sample_freq/sf_down):size(A{1},1);
  A{1} = A{1}(idx_use,:);
  Depth_speed{1} = Depth_speed{1}(idx_use,:);
end

% Animal flat/ascend/descend (fad) definitions and corresponding indices 
% (in number, rather than bool).
flat_thrs = prctile(Depth_speed{1}, 70);
for i_set = 1:2
  id_fad{i_set}{1} = find(abs(Depth_speed{i_set}) <= flat_thrs); % Flat.
  id_fad{i_set}{2} = find(Depth_speed{i_set} < -flat_thrs); % Ascend.
  id_fad{i_set}{3} = find(Depth_speed{i_set} > flat_thrs); % Descend.
end

% Record Acc_2 flat/ascend/descend status.
fad_status = zeros(size(Depth_2));
for i_fad = 1:3
  fad_status(id_fad{2}{i_fad}) = i_fad;
end


%% Cluster A{1} to tight up the distribution, testing, not used.
tight_up_a1 = false;
if tight_up_a1
  epsilon = 0.1;
  minpts = 20;
  for i_fad = 1:3
    i_clus{i_fad} = dbscan(A{1}(id_fad{1}{i_fad},:), epsilon, minpts);
  end
  
  if plot_on
    %% Plot clusters.
    figure
    hold on
    for i_fad = 1:3
      i_use = id_fad{1}{i_fad}(i_clus{i_fad} == 1);
      i_notuse = id_fad{1}{i_fad}(i_clus{i_fad} ~= 1);
      plot3(A{1}(i_use, 1), A{1}(i_use, 2), A{1}(i_use, 3), '.')
      plot3(A{1}(i_notuse, 1), A{1}(i_notuse, 2), A{1}(i_notuse, 3), 'm.')

    end
    grid on
    axis equal
    title('Denser A1 cluster')
  end
  % Update id_fad.
  for i_fad = 1:3
    id_fad{1}{i_fad} = id_fad{1}{i_fad}(i_clus{i_fad} == 1);
  end
end

%% Calculate the nearest distance from each point of A2 to the A1 cluster.
% And find the indices of points of A2 that are close to A1.
K = 30;
dist_thrs = 0.1;
id_close = [];
%disp('Finding inliers between sections...')
for i_fad = 1:3 % fad - Flat, Ascend, Descend.
  [dmvar01, dist_mat] = knnsearch(...
    A{1}(id_fad{1}{i_fad}, :),...
    A{2}(id_fad{2}{i_fad}, :),...
    'K', K, 'Distance', 'euclidean');
  dist = median(dist_mat, 2);
  id_close = [id_close;
              id_fad{2}{i_fad}(dist < dist_thrs)];
end
%disp('Inliers found.')
inliers = false(length(A{2}),1);
inliers(id_close) = true;


if plot_on
  %% Plot data.
  for i_set = 1:2
    idx_plot{i_set} = false(length(A{i_set}), 1);
    idx_plot{i_set}(1:plot_interval:length(A{i_set})) = true;
  end

  figure
  i_set = 1;
  %plot3(A{i_set}(idx_plot{i_set}, 1),A{i_set}(idx_plot{i_set}, 2),...
  %  A{i_set}(idx_plot{i_set}, 3),'.')
  plot3(A{i_set}(:, 1),A{i_set}(:, 2),...
    A{i_set}(:, 3),'.')
  hold on
  
  i_set = 2;
  idx_plot_in = idx_plot{i_set} & inliers;
  plot3(A{i_set}(idx_plot_in, 1),A{i_set}(idx_plot_in, 2),...
    A{i_set}(idx_plot_in, 3),'.')
  
  idx_plot_out = idx_plot{i_set} & ~inliers;
  plot3(A{i_set}(idx_plot_out, 1),A{i_set}(idx_plot_out, 2),...
    A{i_set}(idx_plot_out, 3),'.')
  hold off
  
  grid on
  axis equal
  xlabel('x-surge')
  ylabel('y-sway')
  zlabel('z-heave')
  title('A1 vs. A2 in/out-liers')
  legend('A1', 'A2-in', 'A2-out')

end

