function [sinuousity,range,speed]= getsinuousity6(x,y,Fs)
% calculate sinuousity speed and range across x and y line segment
% assumes Fs of 10 Hz in x and y so output will have to be resampled
% returns one value which is sinuosity,range and speed determined across
% whole block

% based on getsinuosity5   17/5/24


% calculate the duration of the line segment
duration= length(x)/Fs;         % duration of block in seconds

% Calculate the length of the line segment
segment_length = sqrt(diff(x).^2 + diff(y).^2);         % separate value for each increment
total_length=nansum(segment_length);

% calculate mean speed across block
speed=total_length/duration;

% Calculate the curvature of the line segment for every block
span=1;   % how many points to compare across
curvature = zeros(size(x));
for i = span+1:length(x)-span
    % Calculate the distance between consecutive points
    d1 = sqrt((x(i) - x(i-span))^2 + (y(i) - y(i-span))^2);
    d2 = sqrt((x(i+span) - x(i))^2 + (y(i+span) - y(i))^2);
    
    % Calculate the angle between the vectors formed by the adjacent points
    dot_product = (x(i+span) - x(i))*(x(i) - x(i-span)) + (y(i+span) - y(i))*(y(i) - y(i-span));
    angle = acos(dot_product / (d1 * d2));
    
    % Calculate the curvature
    curvature(i) = 2 * sin(angle) / segment_length(i);
end
sinuousity=mean(curvature);   % mean for this block

% Calculate the eccentricity and range of the line segment for every block

[eccentricity,range] = get_path_eccentricity(x,y);

