function inlier_density_plt(TagData)
%INLIER_DENSITY_PLT Helper function to plot inlier densities for auto orient

T_m = (1:length(TagData.accelTag))/TagData.sampleFreq/60;
in_density_forward = TagData.in_density_forward;
in_density_reverse = TagData.in_density_reverse;
density_thrs = TagData.autoOrientPars.pars.density_thrs;
idx_breaks = TagData.autoOrientPars.idx_breaks;
sec_idx_list = TagData.autoOrientPars.sec_idx_list;

figure
axd(1) = subplot(2,1,1);
hold on
plot(T_m, in_density_forward(:,1), 'lineWidth',1)
plot(T_m, in_density_forward(:,2), 'lineWidth',2)
plot(T_m([1,end]),[density_thrs, density_thrs], 'm:', 'lineWidth',1.5)

% Plot one tag move time and sub section break first to easy legend
% adding.
i = 1;
i_slide = sec_idx_list(i);
plot([T_m(i_slide),T_m(i_slide)], [0,1],'k-.', 'lineWidth',2)
i_start = idx_breaks(i);
plot([T_m(i_start),T_m(i_start)], [0,1],'r:')

for i = 2:length(sec_idx_list)
    i_slide = sec_idx_list(i);
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
for i = 1:length(sec_idx_list)
    i_slide = sec_idx_list(i);
    plot([T_m(i_slide),T_m(i_slide)], [0,1],'k-.', 'lineWidth',2)
end

%grid on
xlabel('time (m)')
ylabel('inlier density')
title('Checking backward')
linkaxes(axd, 'x')

end