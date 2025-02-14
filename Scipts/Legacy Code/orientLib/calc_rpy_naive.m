function [roll, pitch, yaw, A_g, A_spec] = ...
  calc_rpy_naive(A, M, sample_freq, wd_naive)
% Estimate roll, pitch and yaw of the animal assuming accelerometer
% reading is all about gravity. The input sensor data should be aligned
% with the animal's body frame: x-surge points forward, y-sway points left,
% and z-sway points up. If input is not aligned with the animal's body,
% then this function estimates the rpy of the tag itself rather than the 
% animal.
%
% INPUT:
% A             n-by-3 acceleromter reading [unit does not matter]. 
%               Axes are aligned with the animal's body frame as defined 
%               above.
% M             n-by-3 *calibrated* magnetometer reading [unit does not matter]. 
%               Axes are aligned with the animal's body frame as defined 
%               above.
% sample_freq   1-by-1 scalar, the sampling frequency of the sensor.
% wd_naive      1-by-1 scalar, window size of the (low-pass) moving
%               average filter, used in the naive method for identifying 
%               portion of the signal associated static and dynamic
%               orientation, fitering out linear accelrations. 
%               So the filtered signal is getting closer to the assumption 
%               accelerometer reading is all about gravity. [counts]
%               To start with try: "wd_naive = 0.5 * sample_freq;" 
%
%
% OUTPUT:
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
% A_g           n-by-3 accleration measure (static + dynamic) associated 
%               with gravity. [g]
% A_spec        n-by-3 accleration measure associated with (dynamic) linear
%               accleration of the animal. [g]
%
% ====================
% Ding Zhang
% zhding@umich.edu 
% Updated: 04/06/2021
% ====================

disp('Calculating roll, pitch and yaw using naive approach...')
A = movmean(A, ceil(0.1*sample_freq));
M = movmean(M, ceil(0.1*sample_freq));

% Normalize to 1 g.
% A = A./sqrt(sum(A.^2, 2));
A = A/mean(sqrt(sum(A.^2, 2)));

% Acceleration measure associated with gravity, determined by static and 
% dynamic orientation of the animal.
A_g = movmean(A, wd_naive);

A_spec = A - A_g;

roll = atan2d(A_g(:,2), A_g(:,3));
% Note again that positive rotation around y-axis in our axes definatio
% n is actually negative pitch. 
pitch = real(asind(A_g(:,1)));
yaw = zeros(size(roll));

% Origin position (0 roll and pitch) gravity reading.
%A_g = [0, 0, 1];
for i = 1:length(yaw)
  % Change of coordinate to find the equivalent magnetic measurement as 
  % if the tag is flat (i.e. zero roll and pitch).
  %R = vrrotvec2mat(vrrotvec(A_g, A(i,:))); % Not stable.
  R = rot_xd(roll(i)) * rot_yd(-pitch(i)); % Note again pitch is negated.
  m = M(i,:)*R;
%   try
  yaw(i) = -atan2d(m(2),m(1));
%   catch
%     disp('here')
%   end
end

disp('Roll, pitch and yaw calculation done.')

%%
% % Testing change of coordinates.
% A3 = zeros(size(A));
% A5 = zeros(size(A));
% A_g = [0, 0, 1];
% for i = 1:length(roll)
%   R3 = rot_xd(roll(i)) * rot_yd(-pitch(i)); 
%   A3(i,:) = A(i,:)*R3; % A3 -> A_g
%   R5 = rot_yd(pitch(i)) * rot_xd(-roll(i));
%   A5(i,:) = A_g*R5; % A5 -> A
% end

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