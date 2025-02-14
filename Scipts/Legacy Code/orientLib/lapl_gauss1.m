function [x, y] = lapl_gauss1(sig, fs)
%LAPL_GAUSS1 Computes laplacian of gaussian to produce 1d kernel
%   Also uses sampling frequency for point interval scaling

% Generate input vec
bnd = ceil(4*sig);
x = (-bnd:(1/fs):bnd)';

% Compute LoG kernel
% NOTE: SIGMA EXPONENT IN DENOMINATOR IS 3, BUT ORIGINAL DERIVATIVE IS 5, THIS
% IS BECAUSE KERNEL NEEDS TO BE SCALED BY SIGMA^2 TO ENSURE SCALE IS NORMALIZED
% WHEN RUN THROUGH CONVOLUTION
y = ((x.^2 - sig.^2)/((sig.^3)*sqrt(2*pi))).*exp(-x.^2/(2*sig.^2));

end