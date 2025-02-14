function Ry = rot_yd(beta)
% Rotation matrix around y axis. Input 'beta' angle in degree.
% Right side multiplication: "v2 = v1*Ry" with 1-by-3 v1 and v2.
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 05/21/2021
% ====================
Ry = [cosd(beta),   0, sind(beta);
      0,            1,          0;
      -sind(beta),  0, cosd(beta)]';
end

