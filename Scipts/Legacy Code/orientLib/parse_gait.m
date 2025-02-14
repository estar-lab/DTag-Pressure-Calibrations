function [frequency_filt, amplitude_filt, period_filt, bool_active,...
  frequency_at_peaks, amplitude_at_peaks, period_at_peaks,...
  positive_peaks, negative_peaks] =...
  parse_gait(signal, sample_freq, min_amplitude, max_period,...
  plot_on, signal_name, signal_unit)
% Parse the gait of the animal using the input signal (e.g. pitch
% of dynamic pose) by locating positive and negative peaks in signal, 
% calculate gait parameters (stroke frequncy, amplitude, period) 
% and identify active swimming data segments. 
% The data segments that do not contain peaks would be considered inactive
% (e.g. gliding).
%
% INPUT:
% signal          n-by-1 input signal, e.g. pitch of dynamic pose or simply
%                 pitch, for example.
% sample_freq     1-by-1 scalar of sample frequncy of the input signal. [Hz]
% min_amplitude   1-by-1 scalar indicates the minimum amplitude
%                 (peak-to-peak) for a qualified peak in signal.
% max_period      1-by-1 scalar indicates the maximum period for a qualified 
%                 peak in signal.
% plot_on         1-by-1 boolean variable to trun on or off plot.
% signal_name     Chars specify the name of signal, for plot only.
% signal_unit     Chars specify the unit of signal, for plot only.
%
% OUTPUT:
% frequncy_filt   n-by-1 stroke frequncy of the animal.
% amplitude_filt  n-by-1 stroke amplitude (peak-to-peak) of the animal.
% period_filt     n-by-1 stroke period of the animal.
% bool_active     n-by-1 boolean indicates whether a point is associated
%                 with active swimming. 
%
% frequncy_at_peaks   n-by-1 stroke frequncy calculated at each peak, 
%                     values are NaN for points that are not peak.
% amplitude_at_peaks  n-by-1 stroke amplitude calculated at each peak, 
%                     values are NaN for points that are not peak.  
% period_at_peaks     n-by-1 stroke period calculated at each peak, 
%                     values are NaN for points that are not peak.
%
% positive_peaks  n-ny-1 values for the detected positive peaks of the
%                 signal, values are NaN for points that are not peak.
% negative_peaks  n-ny-1 values for the detected negative peaks of the
%                 signal, values are NaN for points that are not peak.
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 04/11/2021
% ====================

if nargin == 2
  min_amplitude = 6; % [degree]
  max_period = 3; % [s].
  plot_on = false;
  signal_name = 'pitch_{dp}';
  signal_unit = 'deg';
elseif nargin == 4
  plot_on = false;
  signal_name = 'pitch_{dp}';
  signal_unit = 'deg';
elseif nargin == 5
  signal_name = 'pitch_{dp}';
  signal_unit = 'deg';
end


%mpp = min_amplitude; % Min Peak Prominence. [degree]
mpw = round(max_period/2*sample_freq); % Max Peak Width.
[pks_p, locs_p, w_p, p_p] = findpeaks(signal,...
                                      'MinPeakProminence',min_amplitude,...
                                      'MaxPeakWidth', mpw);
[pks_n, locs_n, w_n, p_n] = findpeaks(-signal,...
                                      'MinPeakProminence',min_amplitude,...
                                      'MaxPeakWidth', mpw);
% Active swim detection.
bool_active = false(size(signal));
bool_active(abs(signal) > min_amplitude/2) = true;

for i = 1:length(locs_p)
  i_s = max(1, round(locs_p(i) - w_p(i)));
  i_e = min(length(signal), round(locs_p(i) + w_p(i)));
  bool_active(i_s:i_e) = true;
end
for i = 1:length(locs_n)
  i_s = max(1, round(locs_n(i) - w_n(i)));
  i_e = min(length(signal), round(locs_n(i) + w_n(i)));
  bool_active(i_s:i_e) = true;
end

% Positive and negative peaks.
positive_peaks = nan(size(signal));
positive_peaks(locs_p) = pks_p;
negative_peaks = nan(size(signal));
negative_peaks(locs_n) = -pks_n;

% (Peak-to-Peak) Amplitude calculation.                                    
amplitude_at_peaks = nan(size(signal));
amplitude_at_peaks(locs_p) = p_p;
amplitude_at_peaks(locs_n) = p_n;

%amplitude_filt = movmean(amplitude_at_peaks, mpw,'omitnan');
amplitude_filt = movmedian(amplitude_at_peaks, mpw,'omitnan');
amplitude_filt(~bool_active) = 0;
amplitude_filt(isnan(amplitude_filt)) = 0;

% Period calculation.
period_at_peaks = nan(size(signal));
period_at_peaks(locs_p) = 2*w_p/sample_freq;
period_at_peaks(locs_n) = 2*w_n/sample_freq;

%period_filt = movmean(period_at_peaks, mpw,'omitnan');
period_filt = movmedian(period_at_peaks, mpw,'omitnan');
period_filt (~bool_active) = 0;
period_filt (isnan(period_filt)) = 0;

% Frequncy calculation.
frequency_at_peaks = 1./period_at_peaks;

frequency_filt = movmean(frequency_at_peaks, mpw,'omitnan');
frequency_filt(~bool_active) = 0;
frequency_filt(isnan(frequency_filt)) = 0;


%% Plot peaks.
if plot_on
  time = (1:length(signal))/sample_freq;
  signal_name = convertStringsToChars(signal_name);
  signal_unit = convertStringsToChars(signal_unit);

  figure
  m_pk = 3;
  ax_pk(1) = subplot(m_pk,1,1);
  plot(time, signal)
  hold on
  plot(time(locs_p), pks_p, '^')
  plot(time(locs_n), -pks_n, 'v')
  plot(time(~bool_active), signal(~bool_active), '.')
  grid on
  xlabel('time (s)')
  ylabel([signal_name, ' (', signal_unit, ')'])
  legend(signal_name, 'positive peaks', 'negative peaks', 'inactive')

  ax_pk(2) = subplot(m_pk,1,2);
  plot(time, amplitude_at_peaks, 'x')
  hold on
  plot(time, amplitude_filt, '-')
  grid on
  xlabel('time (s)')
  ylabel(['amplitude (', signal_unit, ')'])
  legend('peaks', 'filtered')

  ax_pk(3) = subplot(m_pk,1,3);
  plot(time, frequency_at_peaks, 'o')
  hold on
  plot(time, frequency_filt, '-.')
  grid on
  xlabel('time (s)')
  ylabel('frequncy (Hz)')
  legend('peaks', 'filtered')

  % ax_pk(4) = subplot(m_pk,1,4);
  % plot(time, period_stack, 'x')
  % hold on
  % plot(time, period_filt, '-')
  % 
  % grid on
  % xlabel('time (s)')
  % ylabel('period (s)')
  % 
  linkaxes(ax_pk, 'x')
  
end
