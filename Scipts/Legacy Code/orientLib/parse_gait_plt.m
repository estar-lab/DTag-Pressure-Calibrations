function parse_gait_plt(TagData)
%PARSE_GAIT_PLT Helper function to plot results of parse_gait

sample_freq = TagData.sampleFreq;
signal = TagData.relPitch;
signal_name = TagData.autoOrientPars.signal_name;
signal_unit = TagData.autoOrientPars.signal_unit;
positive_peaks = TagData.upFluke;
negative_peaks = TagData.dnFluke;
locs_p = find(~isnan(positive_peaks));
locs_n = find(~isnan(negative_peaks));
pks_p = positive_peaks(~isnan(positive_peaks));
pks_n = -negative_peaks(~isnan(negative_peaks));
bool_active = ~TagData.inactive;
amplitude_at_peaks = TagData.atpeakAmp;
amplitude_filt = TagData.filtAmp;
frequency_at_peaks = TagData.atpeakFreq;
frequency_filt = TagData.filtFreq;

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
ylabel(['fluke amp. (', signal_unit, ')'])
legend('peaks', 'filtered')

ax_pk(3) = subplot(m_pk,1,3);
plot(time, frequency_at_peaks, 'o')
hold on
plot(time, frequency_filt, '-.')
grid on
xlabel('time (s)')
ylabel('fluke freq. (Hz)')
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