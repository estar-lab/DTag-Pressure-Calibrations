function auto_orient_res_plt(TagData)
%AUTO_ORIENT_RES_PLT Plots summary results from Auto Orientation estimator

% Extract parameters
sample_freq = TagData.sampleFreq;
A_org = TagData.accelInterm;
Depth = TagData.depth;
roll_niv = TagData.roll;
pitch_niv = TagData.pitch;
yaw_niv = TagData.head;
signal = TagData.relPitch;
positive_peaks = TagData.upFluke;
negative_peaks = TagData.dnFluke;
bool_active = ~TagData.inactive;
amplitude_at_peaks = TagData.atpeakAmp;
amplitude_filt = TagData.filtAmp;
frequency_at_peaks = TagData.atpeakFreq;
frequency_filt = TagData.filtFreq;


% Plot
time = (1:length(A_org))'/(sample_freq);
time_unit = 's';

figure
m_pk = 5;

lw = 1.2;
lw2 = 1.2;

font_size = 12;

i_p = 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
hold on
plot(time, Depth, 'k-', 'lineWidth', lw)
grid on
ylabel('Depth (m)', 'FontSize', font_size)


i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
hold on
plot(time, roll_niv, 'k-', 'lineWidth', lw)
plot(time, pitch_niv, 'b-', 'lineWidth', lw)
plot(time, yaw_niv, 'r-', 'lineWidth', lw)
grid on
ylabel('Orientation (deg)', 'FontSize', font_size)
% legend('roll', 'pitch', 'yaw',...
%   'FontSize', font_size, 'Orientation','horizontal')
legend('roll', 'pitch', 'yaw', 'FontSize', font_size)
% legend('boxoff')


i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
plot(time, signal, 'k-', 'lineWidth', lw2)
hold on
plot(time, positive_peaks, '^')
plot(time, negative_peaks, 'v')
plot(time(~bool_active), signal(~bool_active), 'r.')

grid on
ylabel('Pitch_{dp} (deg)', 'FontSize', font_size)
% legend('signal', 'positive peak', 'negative peak', 'inactive',...
%   'FontSize', font_size, 'Orientation','horizontal')
legend('signal', 'positive peak', 'negative peak', 'inactive', 'FontSize', font_size)
% legend('boxoff')


i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
plot(time, frequency_at_peaks, 'bo', 'lineWidth', lw2)
hold on
plot(time, frequency_filt, 'k-', 'lineWidth', lw2)
grid on
ylabel('Fluke Freq. (Hz)', 'FontSize', font_size)
% legend('at-peak', 'filtered', 'FontSize', font_size, ...
%   'Orientation','horizontal')
legend('at-peak', 'filtered', 'FontSize', font_size)
% legend('boxoff')


i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
plot(time, amplitude_at_peaks, 'bx', 'lineWidth', lw2)
hold on
plot(time, amplitude_filt, 'k-', 'lineWidth', lw2)
grid on
xlabel(['Time (', time_unit, ')'], 'FontSize', font_size)
ylabel('Fluke Amp. (deg)', 'FontSize', font_size)
% legend('at-peak', 'filtered', 'FontSize', font_size, ...
%   'Orientation','horizontal')
legend('at-peak', 'filtered', 'FontSize', font_size)
% legend('boxoff')

linkaxes(ax_pk, 'x')

end