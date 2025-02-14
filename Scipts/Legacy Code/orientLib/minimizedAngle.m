function angle = minimizedAngle(angle)
% Wraps angle to (-PI,PI)

angle = mod(angle + pi, 2*pi) - pi;

end