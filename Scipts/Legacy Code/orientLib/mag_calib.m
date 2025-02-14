function [M, idx_cal] = mag_calib(M_raw, A, Depth, idx_cal)
% Magnetometer calibration function. First tries to find section of data
% that the tag is under water for mag offset calibration. If few data points
% are found, shows the data to user to select section of data that the tag 
% is away from metal or concrete. Those data will be used to find the offset
% of the magnetometer reading. By fitting a ball to the selected section of
% data, find the center of the ball, and subtract center so the ball is 
% centered around zero.
%
% Input:
% M_raw    Magnetometer reading, m-by-3, going to be calibrated.
%
% A        Accelerometer reading, m-by-3, for user view only.
%
% Depth    Pressure/depth reading, m-by-1, for finding underwater section or
%          for user view.
%
% idx_cal (Optional) The indices for the section of data that is used for
%          magnetometer calibration. 
%          If not provided: the underwater (i.e. abs(depth) > 0.2 meters)
%          section of data will be used for magnetometer calibration.
%          If empty (i.e. idx_cal = []) array is provided: the user will be
%          asked to key it in during the run of the function.
%          If provided (e.g. idx_cal = 1000:50000), then that section will
%          be used to calibrate magnetometer.
%
% Output:
% M     Calibrated magnetometer data, m-by-3.
%
% idx_cal  The indices input by user or found by program for mag calibration.
%
% =====================
% Ding Zhang
% zhding@umich.edu
% Updated: May 30, 2020
% =====================

assert(length(M_raw) == length(A) && length(A) == length(Depth), ...
  'Input data channels dimension mismatch.')

% When plot, do not plot every points.
idx_plot = 1:10:length(A);


if nargin < 4
  % 'idx_cal' is not provided. Use 'Depth' to decide section of data for
  % calibration.
  if mean(Depth) < 0
    idx_cal = find(Depth<-0.2);
  else
    idx_cal = find(Depth>0.2);
  end
end

% If there is not enough calibration points, let the user know.
n_rep = 0;
while length(idx_cal) < 100
  n_rep = n_rep + 1;
  if n_rep > 10
    error('Tried too many times, gyro calib program terminates. Have a good time.')
  end
  disp('Not enough calibration points, length(idx_cal)<100, more is better.')
  disp('With the data plot, specify indices for the section of data to use for ')
  disp('calibration.')
  idx_cal = get_user_input_for_calibration(A, M_raw, Depth, idx_plot);
end  
% Find magnetometer offset. 
M_os_xyz = fit_ball_v2(M_raw(idx_cal,:));

% Shift magnetometer data.
M = M_raw - ones(length(M_raw),1)*M_os_xyz';

% Show calibrated data to user if user input was required.
if n_rep > 0
  figure
  idx_plot_cal = intersect(idx_plot,idx_cal);
  plot3(M(idx_plot_cal,1), M(idx_plot_cal,2), M(idx_plot_cal,3),'.')
  xlabel('x')
  ylabel('y')
  zlabel('z')
  grid on
  axis equal
  title({'Section of magnetometer data after calibration.';
         'Should look like a ball centered around origin of the axes.';
         'Otherwise the oriantation estimation may not be accurate.'})
end
disp('Magnetometer calibration done.')
end


function [center,r] = fit_ball_v2(data)
%This is a solution based on linearizing the sphere equation.
%The solution is not least-squares in its truest form, however this result
%is good enough for most purposes. This can be an input for the initual
%guess while using the Gauss-Newton or the Levenberg?Marquardt algorithm 
%
% The data has to be in 3 columns and at least 4 rows, first column with Xs, 
% 2nd column with Ys and 3rd columns with Zs of the sphere data. 
% The output 
% r = radius
% a = X coordinate of the center
% b = Y coordinate of the center
% c = Z coordinate of teh center
% center = [a;b;c]
%
% usage: fit_ball_v2(data) % where data is a mx3 data and m>=4
% Updated by Ding Zhang Dec 31, 2019.

xx = data(:,1);
yy = data(:,2);
zz = data(:,3);
AA = [-2*xx, -2*yy , -2*zz , ones(size(xx))];
BB = [-(xx.^2+yy.^2+zz.^2)];
YY = mldivide(AA,BB); %Trying to solve AA*YY = BB
a = YY(1);
b = YY(2);
c = YY(3);
D = YY(4); % D^2 = a^2 + b^2 + c^2 -r^2(where a,b,c are centers)
r = sqrt((a^2+b^2+c^2)-D);
center = [a;b;c];
end


function idx_cal = get_user_input_for_calibration(A, M_raw, Depth, idx_plot)
  % Show data, wait for user's input for indices of the section of data
  % for calibration.
    figure
    axx(1) = subplot(3,1,1);
    plot(idx_plot, A(idx_plot,:))
    xlabel('index of data point')
    ylabel('acc (g)')
    grid on
    title({'This is for magnetometer calibration. Find the start and end indices ';
           'of a long section of data that the tag is away from metal or concrete,';
           'then input the indices into the command window.';
           'E.g. input 100:6000 if the start index is 100 and end index is 6000.'})

    axx(2) = subplot(3,1,2);
    plot(idx_plot, M_raw(idx_plot,:))
    xlabel('index of data point')
    ylabel('mag (uT)')
    grid on

    axx(3) = subplot(3,1,3);
    plot(idx_plot, Depth(idx_plot), 'x')
    xlabel('index of data point')
    ylabel('depth (m)')
    grid on
    linkaxes(axx,'x')

    disp('Check the plot and key in the indices for magnetometer calibration')
    prompt = '(e.g. 100:6000): ';
    idx_cal = input(prompt);
end
