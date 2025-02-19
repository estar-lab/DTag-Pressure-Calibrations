clear; close all;

select_data_dir();
rec_names = {'345PressCal_11022025';
             '338PressCal_11022025';
             '344PressCal_11022025';
             '340PressCal_12022025';
             '339PressCal_12022025';
             '348PressCal_12022025'};

for i = 1:length(rec_names)
    % Add DTag library to path
    addpath(genpath('Legacy Code\dtagLib\'), 'Legacy Code\orientLib\');
    
    % Set main data directory
    if ~exist('data_dir.mat', 'file')
        disp('No data directory file found: Run previous section')
    end
    dir_st = load('data_dir.mat');
    data_dir = dir_st.dirname;
    rec_name = rec_names{i}; % Blue whale
    dir_name = sprintf('%s%s', data_dir, rec_name);
    addpath(genpath(dir_name));
    df = 10;
    lpf = [0, 0, 0];
    filt_pars = ([]);
    filt_pars.lpf = lpf;
    
    est_mthd = 'auto';
    
    TagData.RawVolt = dtag3_test(dir_name, rec_name, df);
    M1 = cell2mat(TagData.RawVolt.x(1:3)');
    P = TagData.RawVolt.x{10};
    plot(P(1:2:end)); ylabel('P Volt'); xlabel('Time (min)');
    cont = 0;
    while ~cont
        [idx_offset, dummy_val] = getpts(gca);
        P(1:round(2*idx_offset)) = dummy_val;
        plot(P(1:2:end)); ylabel('P Volt'); xlabel('Time (min)');
        cont = input('continue? ');
        
    end
    close gcf
    time_rv = (0:(length(M1) - 1))/(TagData.RawVolt.fs(1)*60);
    
    press{i} = P;
    t{i} = time_rv;

end

%% 

% Example input signals (replace with your actual signals)
signals = press;  % Cell array of signals
n = length(signals);  % Number of signals

% Parameters for alignment
alignment_point = 1;  % The reference signal index (for example, the first one)

% Pre-allocate an array to store the aligned signals
aligned_signals = cell(n,1);

% Select the reference signal to align others with
ref_signal = signals{alignment_point};

% Create a new figure for visualization
figure;

% Zero Padding
for i = 1:n
    lens(i) = length(signals{i});
end
max_len = max(lens);
min_len = min(lens);

for i = 1:n
    signal = signals{i};
    signal(end+1:max_len+1) = signal(end)*ones(max_len+1-length(signal),1);
    signals{i} = signal;
end

for i = 1:n
    rec_names{i}(end-8) = ' ';
end
% Loop through the signals and align them with the reference signal

for i = 1:n
    % Compute the cross-correlation between the reference signal and the current signal
    [correlation, lag] = xcorr(ref_signal, signals{i});
    
    % Find the lag with the maximum correlation
    [~, max_lag_idx] = max(abs(correlation));
    lag_value = lag(max_lag_idx);  % This is the time shift (lag) to align the signals
    
    % Align the current signal by shifting it according to the lag
    aligned_signals{i} = circshift(signals{i}, lag_value);
    
    % Find movmean
    median_signals{i} = movmedian(aligned_signals{i}, 1e3);
    
    % Plot the aligned signals for visualization
    
    
end



depths = 1e2*(9:-1:0);


for i = 1:n
    cont = 0;
    while ~cont
        plot(signals{i})
        title(rec_names{i})
        [idxs,~] = getpts(gcf);
        raw_vals = signals{i}(round(idxs));
        if ~isequal(length(raw_vals), length(depths))
            fprintf('Number of points is not the same as the number of depths\n')
            close gcf
        else
            cont = 1;
        end
    end
    params = polyfit(raw_vals, depths, 1);
    corr_signals{i} = params(1)*aligned_signals{i}+params(2);
    Params{i} = params;
end
T = tiledlayout(3,1, "Padding","compact", "TileSpacing","compact");
for i = 1:n    
% Plot everything
    ax1 = nexttile(1);
    hold on
    if all(rec_names{i} == '339PressCal 12022025') || all(rec_names{i} == '348PressCal 12022025')
        plot(signals{i}, LineStyle="--", LineWidth=1); 
    else
        plot(signals{i});
    end
    hold off
    
   
    


    ax2 = nexttile(2);
    
    hold on
    if all(rec_names{i} == '339PressCal 12022025') || all(rec_names{i} == '348PressCal 12022025')
        plot(aligned_signals{i}, LineStyle="--", LineWidth=1)
    else
        plot(aligned_signals{i})
    end
    hold off

    

    ax3 = nexttile(3);
    hold on
    if all(rec_names{i} == '339PressCal 12022025') || all(rec_names{i} == '348PressCal 12022025')
        plot(corr_signals{i}, LineStyle="--", LineWidth=1)
    else
       plot(corr_signals{i});
    end
    hold off
    
    
    

end
linkaxes([ax1, ax2, ax3], 'x')
linkaxes([ax1, ax2], 'y')

legend(ax2, rec_names)
ylabel(ax2, 'P (V)')
title(ax2, 'Aligned Signals')

legend(ax1, rec_names)
ylabel(ax1, 'P (V)')
title(ax1, 'Unaligned Signals')

legend(ax3, rec_names)
ylabel(ax3, 'Depth (m)')
title(ax3, 'Corrected Depths')

% Optional
% saveas(gcf, ['Plots\PressCal_' rec_names{1}(end-7:end) '.png'])

delete('data_dir.mat')