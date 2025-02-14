function TagData = preorient_filter(TagData, lpf)
%PREORIENT_FILTER Apply low-pass filter to data before orientation calculations
%   Applies Savitsky-Golay filtering to data. Processes accelerometer,
%   magnetometer, and pressure data individually. Each component of lpf
%   represents the filter window used for each, format of [accel., mag., depth],
%   windows in seconds (recommended <1s).

fs = TagData.sampleFreq;

lpf = lpf/2; % Cut in half to correctly compute half-windows for filt

% Process accelerometer filtering
if lpf(1)
    whf = round(lpf(1)*fs); % Compute half-window
    TagData.accelTagOrig = savGol(TagData.accelTagOrig, whf, whf, 2);
    TagData.accelTag = savGol(TagData.accelTag, whf, whf, 2);
end

% Process magnetometer filtering
if lpf(2)
    whf = round(lpf(2)*fs); % Compute half-window
    TagData.magTagOrig = savGol(TagData.magTagOrig, whf, whf, 2);
    TagData.magTag = savGol(TagData.magTag, whf, whf, 2);
end

% Process depth filtering
if lpf(3)
    whf = round(lpf(3)*fs); % Compute half-window
    TagData.depthOrig = savGol(TagData.depthOrig, whf, whf, 2);
    TagData.depth = savGol(TagData.depth, whf, whf, 2);
end

end