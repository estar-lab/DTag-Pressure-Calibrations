function [section_idx_list, in_density_forward, in_density_reverse, idx_breaks] =...
  find_tag_slide_times_func(Acc, Depth, sample_freq, sec_dur_min,...
  density_thrs, plot_on)
% Find the indices corresponding to tag sliding/shifting times on the
% animal. This function checks the entire dataset incrementally one 
% *subsection* (length of which defined by "sec_dur_min") at a time. 
% If a slide/shift happend within a subsection, then the data distribution
% of that subsection, starting when the slide/shift happens, will be 
% different from earlier subsections. And the data points after the
% slide/shift time will appear to be outliers with respect to the earlier
% subsections. Built on this observation, we find slide/shift by looking
% for the transition time that *inlier density* drops below a given
% threshold ("density_thrs"). Check the code for more detail.
%
% NOTE: The design of this function assumes that there are at most *two*
% slide/shift instances within *two* subsections. Also, the method is a
% distribution based method, and if two slides/shifts are too close to 
% each other, e.g. within 2 minutes (defined by the variable "wd"), 
% the method tends to not able to detect both slides.
%
% INPUT:
% Acc           n-by-3 tag accelerometer data. [unit does not matter]
% Depth         n-by-1 *calibrated* tag depth data. [m] 
% sample_freq   1-by-1 scalar, sampleing frequency of 'Acc' and 'Depth'.
% sec_dur_min   [Optional] 1-by-1 scalar, defines the duration of each 
%               *subsection* in minutes. Default value is 10 minutes, 
%               the function may not work well if this value is set too low.
% density_thrs  [Optional] 1-by-1 scalar, defines the transition threshold 
%               for inlier density. When density drops below this threshold,
%               the instance is considered a slide/shift. Default value is
%               0.4.
% plot_on       [Optional] 1-by-1 binary flag controls if to show inlier
%               density plot. Default value is true.
%
% OUTPUT:
% section_idx_list    (k+2)-by-1 data point indicies of the k detected
%                     slides/shifts and a "1" at the beginning and a 
%                     "length(Acc)" at the end. So "section_idx_list"
%                     always gives some valid data sections that are
%                     sperated by detected slides/shifts.
% in_density_forward  n-by-2 [inlier density, filtered inlier density] 
%                     corresponds to each time instance. The forward
%                     density will be checked over every *subsection*, 
%                     except the first, of a big section. The first  
%                     subsection (together with the second if available) 
%                     is used as reference to check inlier densities of 
%                     all other subsections, until a slide/shift is found.
%                     Note that, a *subsection* is defined by the function 
%                     to incrementally check the dataset. A (big) section is
%                     defined by slide/shifts. Every slide/shift starts
%                     a new (big) section.
% in_density_reverse  n-by-2 [inlier density, filtered inlier density] 
%                     corresponds to each time instance. The reverse
%                     density will only be checked over every *first* 
%                     *subsection* using the second *subsection* of each 
%                     (big) section.
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 07/15/2020
% 04/05/2022: Added a few lines to flip and filter depth.
% ====================

if nargin == 3
  sec_dur_min = 10;%15
  density_thrs = 0.4;
  plot_on = true;
elseif nargin == 4
  density_thrs = 0.4;
  plot_on = true;
elseif nargin == 5
  plot_on = true;
end  


% Flip Depth if they are positive. In this script, all depth (i.e. z-axis
% position) are assumed negative or zero, so when the animal descends, depth
% decreases.
if nanmean(Depth) > 0
  Depth = -Depth;
end

% Filter depth.
Depth = movmean(Depth, ceil(0.5*sample_freq));


% Defines the section size used in the greedy approach.
sec_dur = round(sec_dur_min*60*sample_freq);

% Window size used to calculate inlier density over time. Given binary 
% in-or-out values of individual data points.
%wd = round(120*sample_freq); % second*sample_freq = num_data_points
wd = round(3*60*sample_freq); % second*sample_freq = num_data_points
% A second window used to filter the inlier density resulted from the
% first window filter.
%wd2 = round(150*sample_freq);
wd2 = round(30*sample_freq);


% To record inlier densities. Two columns [density, filtered density]
in_density_forward = nan(length(Depth), 2);
in_density_reverse = nan(length(Depth), 2);

% Define when two time instances are too close to each other, actions are 
% taken at circumstances when time instances are close.
%near_thrs = round(90*sample_freq);
near_thrs = round(180*sample_freq);

% Stores the detected sliding indices.
idx_slide_detect = [];

% Start index of 1st section.
i_s = 1; 
% End index of 1st section and start of 2nd section.
i_m = i_s + sec_dur;
% End index of 2nd section.
i_e = i_m + sec_dur;

% Store subsection breaking indices.
idx_breaks = [i_m];

disp('Finding tag slide times, may take a few seconds...')
tic
% Run until there is no or very little data left for the 2nd section.
while i_m < length(Acc) - near_thrs
  % Check 2nd section with respect to 1st section.
  % If the 1st section is long (i.e. contains multiple subsections
  % already), then we only use the earlier subsection of the 1st section
  % rather than the entire 1st section, so the data points are not too dense
  % and incremental tag slides can be handled to some extents. 
  i_m_temp = min(i_m, i_s + round(1.5*sec_dur));
  inliers_2 = find_inliers_func(...
    Acc(i_s:i_m_temp,:),Depth(i_s:i_m_temp,:),...
    Acc(i_m:i_e,:), Depth(i_m:i_e,:), sample_freq, []); 
  in_dens_2 = movmean(inliers_2, wd, 'omitnan');
  in_dens_filt_2 = movmean(in_dens_2, wd2, 'omitnan');
  i_slide_2 = find(in_dens_filt_2 < density_thrs, 1, 'first');
  
  % Record density.
  in_density_forward(i_m:i_e,1) = in_dens_2;
  in_density_forward(i_m:i_e,2) = in_dens_filt_2;
  
  
  check_1st_section = false;
  if isempty(i_slide_2) 
    % 2nd section is good.
    if i_m - i_s == sec_dur
      % First section has not been check, check first section.
      check_1st_section = true;
    end
    
  else
    % 2nd section is not good.
    if (i_m - i_s == sec_dur) && (i_slide_2 < near_thrs)
      % 1st section has not been check, and slide point in 2nd section is
      % close to the start of 2nd section, check first section.
      check_1st_section = true;
    else
      % Add i_slide_2 as a tag sliding instance.
      idx_slide = i_m + i_slide_2 - 1;
      idx_slide_detect(end+1,1) = idx_slide;
      i_s = idx_slide + 1;
      i_e = i_s + sec_dur - 1;
    end
  end
  
  % Check 1st section with respect to 2nd section.
  if check_1st_section
    inliers_1 = find_inliers_func(Acc(i_m:i_e,:), Depth(i_m:i_e,:),...
      Acc(i_s:i_m,:), Depth(i_s:i_m,:), sample_freq, []);  
    in_dens_1 = movmean(inliers_1, wd, 'omitnan');
    in_dens_filt_1 = movmean(in_dens_1, wd2, 'omitnan');
    i_slide_1 = find(in_dens_filt_1 < density_thrs, 1, 'last');
    
    % Record density.
    in_density_reverse(i_s:i_m,1) = in_dens_1;
    in_density_reverse(i_s:i_m,2) = in_dens_filt_1;
    
    if isempty(i_slide_1)
      % 1st section is good.
      if ~isempty(i_slide_2)
        % But second section is not good at the beginning.
        % Add i_slide_2 as a tag sliding instance.
        idx_slide = i_m + i_slide_2 - 1;
        idx_slide_detect(end+1,1) = idx_slide;
        i_s = idx_slide + 1;
        i_e = i_s + sec_dur - 1;
      end
      
    else
      % 1st section is not good.
      % Add i_slide_1 as a tag sliding instance.
      idx_slide = i_s + i_slide_1 - 1;
      idx_slide_detect(end+1,1) = idx_slide;
      i_s = idx_slide + 1;
      i_e = i_s + sec_dur - 1;
    end
  end
      
  i_m = i_e + 1;
  i_e = min(i_m + sec_dur, length(Acc));
  idx_breaks(end+1) = i_m;
end
idx_breaks(end) = [];

% Merge close detections and fully define sections.
section_idx_list = [1];
for idx = idx_slide_detect'
  if abs(idx - section_idx_list(end)) < near_thrs
    section_idx_list(end) = round((section_idx_list(end) + idx)/2);
  else
    section_idx_list(end+1) = idx;
  end
end
if abs(length(Acc) - section_idx_list(end)) < near_thrs
  section_idx_list(end) = length(Acc);
else
  section_idx_list(end+1) = length(Acc);
end
section_idx_list(1) = 1;

disp(['Tag slide detection done, number of slides found: ',...
  num2str(length(section_idx_list)-2)])
toc


%% Plot inlier density.
if plot_on
  T_m = (1:length(Acc))/sample_freq/60;
  figure
  axd(1) = subplot(2,1,1);
  hold on
  plot(T_m, in_density_forward(:,1), 'lineWidth',1)
  plot(T_m, in_density_forward(:,2), 'lineWidth',2)
  plot(T_m([1,end]),[density_thrs, density_thrs], 'm:', 'lineWidth',1.5)

  % Plot one tag move time and sub section break first to easy legend
  % adding.
  i = 1;
  i_slide = section_idx_list(i);
  plot([T_m(i_slide),T_m(i_slide)], [0,1],'k-.', 'lineWidth',2)
  i_start = idx_breaks(i);
  plot([T_m(i_start),T_m(i_start)], [0,1],'r:')

  for i = 2:length(section_idx_list)
    i_slide = section_idx_list(i);
    plot([T_m(i_slide),T_m(i_slide)], [0,1],'k-.', 'lineWidth',2)
  end
  for i = 2:length(idx_breaks)
    i_start = idx_breaks(i);
    plot([T_m(i_start),T_m(i_start)], [0,1],'r:')
  end
  %grid on
  xlabel('time (m)')
  ylabel('inlier density')
  title('Checking forward')
  legend('inlier density', 'filtered density', 'threshold',...
    'tag move instances', 'subsection breaks')

  axd(2) = subplot(2,1,2);
  hold on
  plot(T_m, in_density_reverse(:,1), 'lineWidth',1)
  plot(T_m, in_density_reverse(:,2), 'lineWidth',2)
  plot(T_m([1,end]),[density_thrs, density_thrs], 'm:', 'lineWidth',1.5)

  for i = 1:length(idx_breaks)
    i_start = idx_breaks(i);
    plot([T_m(i_start),T_m(i_start)], [0,1],'r:')
  end
  for i = 1:length(section_idx_list)
    i_slide = section_idx_list(i);
    plot([T_m(i_slide),T_m(i_slide)], [0,1],'k-.', 'lineWidth',2)
  end

  %grid on
  xlabel('time (m)')
  ylabel('inlier density')
  title('Checking backward')
  linkaxes(axd, 'x')
end

