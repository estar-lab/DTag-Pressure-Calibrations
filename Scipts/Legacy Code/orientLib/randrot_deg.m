function R = randrot_deg(angle)
% Generate a random rotation matrix that rotates angle degrees in a random
% direction.
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 05/21/2021
% ====================
a_x = [1, 0, 0]; % x-axis, which we used to rotate the fixed degree.
a_rot = randn(1, 3); % The random axis we are rotating around.
% Align a_rot with ax, rotate "angle" degrees around a_x, put a_x back to
% a_rot.
R_rot2x = vrrotvec2mat(vrrotvec(a_rot', a_x'))';
R = R_rot2x*rot_xd(angle)*R_rot2x';
