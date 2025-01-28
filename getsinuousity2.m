function [sinuousity,range,speed]= getsinuousity2(x,y,Fs)
% calculate sinuousity speed and range across x and y line segment
% assumes Fs of 10 Hz so output will have to be resampled

% calculate the duration of the line segment
duration= length(x)/Fs;         % duration in seconds

% Calculate the length of the line segment
segment_length = sqrt(diff(x).^2 + diff(y).^2);

% calculate mean speed across segment
speed=nanmean(segment_length/duration);

% Calculate the curvature of the line segment for every block
curvature = zeros(size(x));
for i = 2:length(x)-1
    % Calculate the distance between consecutive points
    d1 = sqrt((x(i) - x(i-1))^2 + (y(i) - y(i-1))^2);
    d2 = sqrt((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2);
    
    % Calculate the angle between the vectors formed by the adjacent points
    dot_product = (x(i+1) - x(i))*(x(i) - x(i-1)) + (y(i+1) - y(i))*(y(i) - y(i-1));
    angle = acos(dot_product / (d1 * d2));
    
    % Calculate the curvature
    curvature(i) = 2 * sin(angle) / segment_length(i);
end
sinuousity=mean(curvature);   % mean for this block


%% new range version
% Calculate distance of each point from the center
dist2cent = sqrt((x - mean(x)).^2 + (y - mean(y)).^2);
% Calculate total distance from centre along the trajectory
range = sum(dist2cent);
% divide by speed to normalise so that a small range with high speed is a
% lower value than a small range with low speed - small range = exploit

range= range/speed;
range=smooth(range,500);

%% plot if required
ploton=false;
% plot data 
if ploton
    figure
    plot(x,y);
    xlim([-50 50]);
    ylim([-50 50]);
    t=strcat('sinuosity = ',num2str(sinuousity),'   range=', num2str(range)); 
    title(t);
    pause
    close all
end