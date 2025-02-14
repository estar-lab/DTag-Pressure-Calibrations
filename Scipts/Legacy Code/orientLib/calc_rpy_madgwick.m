function [roll, pitch, yaw, A_g, A_spec] = calc_rpy_madgwick(A, G, M, sample_freq)
% Perform Madgwick heading filtering, obtaining orientation estimation and
% acceleration decomposition. Requires 'MadgwickAHRS', 'quaternion_library'
% in the path. 
%
% INPUT:    (Make sure that sensor axis definition are consistent and clear!)
%
% 'A'       n-by-3 Accelerometer reading. [m/s^2 or g] 
%           [unit does not matter for Madgwick, to be normalized]. 
% 'G'       n-by-3 Gyroscope reading, needs to be calibrated with no bias.
%           [rad] (has to be)
% 'M'       n-by-3 Magnetometer reading, needs to be calibrated with no bias.
%           [unit does not matter, to be normalized]
% 'sample_freq'   Sample frequency of sensor.
%
% OUTPUT:
%
% 'A_spec'  Specific acceleration, i.e. the dynamic component of
%           acceleration that is caused by animal's motion directly. Unit
%           is the same as 'A', either m/s^2 or g.
% 'A_g'     Gravity component of the accelerometer reading. (A=A_spec+A_g)
%           Unit is the same as 'A', either m/s^2 or g.
% 'roll'    roll  = euler(1) [rad]
% 'pitch'   pitch = euler(2) [rad]
% 'yaw'     yaw   = euler(3) [rad]
% 'euler'   Describing the tag's pose with-respect-to the "origin" pose,
%           which is defined as [x pointing north, z pointing up, y 
%           pointing west]. ZYX Euler angles where phi--euler(1) is a
%           rotation around X, theta--euler(2) around Y and psi--euler(3)
%           around Z.
% 'rot_mat' Rotation matrices, projects tag frame vector to world frame.
%
% Initial quaternion describing the Earth relative to the sensor. Will
% converge to pose w.r.t. "origin" pose reference within 2 sec.
% "origin" pose is defined as [x pointing north, z pointing up, y pointing
% west]. So 'quaternion_initial' does not have to be changed.
%
%
% Reference:
% [1] Madgwick, S.O., Harrison, A.J. and Vaidyanathan, R., 2011, June. 
% Estimation of IMU and MARG orientation using a gradient descent algorithm. 
% In 2011 IEEE international conference on rehabilitation robotics 
% (pp. 1-7). IEEE.
%
% =======================
% Ding Zhang
% zhding@umich.edu
% Last Updated: 5/05/2020
% =======================

disp('Begin Madgwick orientation filtering...')
quaternion_initial = [1 0 0 0]; 
% Madgwick algorithm gain.
%Beta = 1.6;
%Beta = 0.041;
Beta = 0.1;

AHRS = MadgwickAHRS('SamplePeriod', 1/sample_freq, 'Quaternion',...
    quaternion_initial, 'Beta', Beta);
quaternion = zeros(length(G), 4);
% Main algorithm loop.
for t = 1:length(G)
    AHRS.Update(G(t,:), A(t,:), M(t,:));	% gyroscope units must be radians
    quaternion(t, :) = AHRS.Quaternion;
end
euler = quatern2euler(quaternConj(quaternion)); %
%euler_tag = quatern2euler(quaternion);

roll = rad2deg(euler(:,1));
pitch = -rad2deg(euler(:,2));
yaw = rad2deg(euler(:,3));

% Create array of gravity vectors to remove
g_apx = mode(sqrt(sum(A.^2, 2)));
% Decide unit of A.
if abs(g_apx - 9.8) > abs(g_apx - 1)
  % A is in g.
  g = 1; 
else
  % A is in m/s^2
  g = 9.80297286843;
end
A_g_static = [0, 0, g];
A_g = zeros(size(A));
rot_mat = zeros(3,3,length(A));

for i = 1:length(A)
    A_g(i,:) = A_g_static*quatern2rotMat(quaternConj(quaternion(i,:)));
    rot_mat(:,:,i) = quatern2rotMat(quaternion(i,:));
end

% Subtract gravity array to find specific acceleration
A_spec = A/g_apx*g - A_g;

disp('Orientation filtering done.')


