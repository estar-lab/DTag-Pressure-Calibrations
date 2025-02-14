function X = dtag3_test(recdir, thisDeploy, DF)
%DTAG3_TEST Exports raw dtag3 data

% addpath(genpath(recdir));

if exist([recdir, '/', 'audit'], 'dir') ~= 7
    mkdir([recdir, '/', 'audit'])
end
if exist([recdir, '/', 'cal'], 'dir') ~= 7
    mkdir([recdir, '/', 'cal'])
end
if exist([recdir, '/', 'prh'], 'dir') ~= 7
    mkdir([recdir, '/', 'prh'])
end
if exist([recdir, '/', 'raw'], 'dir') ~= 7
    mkdir([recdir, '/', 'raw'])
end

initPath(recdir)

TagData.deployName = thisDeploy;
TagData.desampFreq = DF;

prefix = TagData.deployName;
X = d3readswv(recdir, prefix, TagData.desampFreq);

end