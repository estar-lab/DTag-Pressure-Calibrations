function pre_orient_correct_plt(TagData, choice)
%PRE_ORIENT_CORRECT_PLT Plots orientation info before corrections are applied

A_org = TagData.accelInterm;
Depth = TagData.depth;
sample_freq = TagData.sampleFreq;
correction_method = TagData.autoOrientPars.pars.correction_method;
sec_idx_list = TagData.autoOrientPars.sec_idx_list;

surf_depth = []; 
plot_interval = [];
if nargin < 2
    choice = 'section'; % 'all_data' or 'section'
end

switch choice
    case 'all_data'
        find_tag_orientation_func(A_org, Depth, sample_freq,...
            surf_depth, plot_interval, correction_method);
    case 'section'
        i_sec = 1;
        i_s = sec_idx_list(i_sec);
        i_e = sec_idx_list(i_sec+1);
        i_use = i_s:i_e;
        find_tag_orientation_func(A_org(i_use,:), Depth(i_use,:),...
        sample_freq, surf_depth, plot_interval, correction_method);
end


end