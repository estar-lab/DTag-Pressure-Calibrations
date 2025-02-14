function select_data_dir()
%SET_DATA_DIR Opens a GUI to help user select the data directory

% Set instructions string
instr = ['Select main tag data directory: Open directory and click ', ...
    '"Select Folder"'];

% Get directory name and print instructions in terminal
fprintf([instr, '\n']);

dirname = uigetdir('...\DTag-Pressure-Calibrations\Data\New Cal Data', instr);

% Break if no folder selected
if ~dirname
    disp('Ended: No directory selected')
    return
end

% Append '\' to close folder path
dirname = [dirname, '/'];

% Set directory in mat file
save('data_dir.mat', 'dirname')

end