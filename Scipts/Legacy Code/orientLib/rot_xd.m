function Rx = rot_xd(alpha)
% Rotation matrix around x axis. Input 'beta' angle in degree.
% Right side multiplication: "v2 = v1*Rx" with 1-by-3 v1 and v2.
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 05/21/2021
% ====================
Rx = [1, 0,           0;
      0, cosd(alpha), -sind(alpha);
      0, sind(alpha), cosd(alpha)]';
end