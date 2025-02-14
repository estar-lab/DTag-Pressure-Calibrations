function Rz = rot_zd(gamma)
% Rotation matrix around z axis. Input 'gamma' angle in degree.
% Right side multiplication: "v2 = v1*Rz" with 1-by-3 v1 and v2.
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 05/21/2021
% ====================
Rz = [cosd(gamma), -sind(gamma), 0;
      sind(gamma), cosd(gamma),  0;
      0,           0,            1]';
end