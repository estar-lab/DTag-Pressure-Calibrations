function [P_dp, Y_dp, V_dp, roll_filt, pitch_filt, yaw_filt] =...
  calc_dynamic_pose(roll, pitch, yaw, wd)
% Estimate roll, pitch and yaw of the animal assuming accelerometer
% reading is all about gravity. The input sensor data should be aligned
% with the animal's body frame: x-surge points forward, y-sway points left,
% and z-sway points up. If input is not aligned with the animal's body,
% then this function estimates the rpy of the tag itself rather than the 
% animal.
%
% INPUT:
% roll          n-by-1 roll angle of the animal in [degree].
%               Roll is defined as the rotation around body forward axis
%               (x-axis). At zero roll, the y-axis has zero gravity
%               reading. Positive roll is body lean to the right.
% pitch         n-by-1 pitch angle of the animal in [degree].
%               Pitch is defined as the angle between body forward axis 
%               (x-axis) and the horizontal plane of the world. 
%               Positive pitch is head up. Note that positive rotation
%               around y-axis in our axes defination is negative pitch. 
% yaw           n-by-1 yaw angle of the animal in [degree]. The direction
%               of magnetic north is yaw 0. West is 90 deg. South is 180.
%               East is -90 deg.
% wd            Window size of the moving average filter. For finding
%               static roll, pitch and yaw.
%
% OUTPUT:
% P_dp          n-by-1 pitch of dynamic pose of the animal in [degree].
% Y_dp          n-by-1 yaw of dynamic pose of the animal in [degree].
% V_dp          n-by-3 3d vector of dynamic pose of the animal.
% roll_filt     n-by-1 low-pass filtered roll angle of the animal in [degree].
% pitch_filt    n-by-1 low-pass filtered pitch angle of the animal in [degree].
% yaw_filt      n-by-1 low-pass filtered yaw angle of the animal in [degree].
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 03/29/2021
% ====================

disp('Calculating dynamic pose...')
wd = max(1, round(wd));
% Note again that positive rotation around y-axis in our axes defination
% is actually negative pitch of the animal.
roll_filt = rad2deg(wrapToPi(movmean(unwrap(deg2rad(roll)), wd)));
pitch_filt = rad2deg(wrapToPi(movmean(unwrap(deg2rad(pitch)), wd)));
%yaw_filt = rad2deg(wrapToPi(movmean(unwrap(deg2rad(yaw)), wd)));
yaw_filt = rad2deg(wrapToPi(movmedian(unwrap(deg2rad(yaw)), wd)));
v0 = [1, 0, 0];
V_dp = zeros(length(roll), 3);
P_dp = zeros(size(pitch));
Y_dp = zeros(size(pitch));

for i = 1:length(yaw)
  R = rot_yd(-pitch(i))*rot_zd(yaw(i))*...
    rot_zd(-yaw_filt(i))*rot_yd(pitch_filt(i))*rot_xd(-roll_filt(i));
  
  %v1 = v0*R;
  V_dp(i, :) = v0*R;
  Y_dp(i) = asind(V_dp(i, 2));
  %v2 = v1*rot_zd(-Y_dp(i));
  P_dp(i) = asind(V_dp(i, 3));
end

disp('Dynamic pose calculation done.')
end


%% Helper rotation functions.
function Rx = rot_xd(alpha)
% Rotation matrix around x axis. Input 'beta' angle in degree.
% Right side multiplication: "v2 = v1*Rx" with 1-by-3 v1 and v2.
Rx = [1, 0,           0;
      0, cosd(alpha), -sind(alpha);
      0, sind(alpha), cosd(alpha)]';
end

function Ry = rot_yd(beta)
% Rotation matrix around y axis. Input 'beta' angle in degree.
% Right side multiplication: "v2 = v1*Ry" with 1-by-3 v1 and v2.
Ry = [cosd(beta),   0, sind(beta);
      0,            1,          0;
      -sind(beta),  0, cosd(beta)]';
end

function Rz = rot_zd(gamma)
% Rotation matrix around z axis. Input 'gamma' angle in degree.
% Right side multiplication: "v2 = v1*Rz" with 1-by-3 v1 and v2.
Rz = [cosd(gamma), -sind(gamma), 0;
      sind(gamma), cosd(gamma),  0;
      0,           0,            1]';
end