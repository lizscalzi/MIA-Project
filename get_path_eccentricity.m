function [eccentricity,range] = get_path_eccentricity(x,y)

% Determine eccentricity and range of path block
% by fitting elipse to path and summing major and minor axes

ploton = false;         % turn on to display results

if ploton
    plot(x, y, 'o-');
    xlabel('x');
    ylabel('y');
    axis equal;
    xlim([-50 50]);
    ylim([-50 50]);
end % ploton

%% Fit an ellipse to the data points
% Prepare data matrix for the fitting
data = [x(:) y(:)];

% Calculate the mean of the data points
mean_data = mean(data);

% Subtract the mean to center the data
centered_data = data - mean_data;

% Perform Singular Value Decomposition (SVD) to find the ellipse axes
[U, S, V] = svd(centered_data, 'econ');

% Semi-major axis length (largest singular value)
a = S(1, 1);

% Semi-minor axis length (smallest singular value)
b = S(2, 2);

% Calculate the eccentricity
eccentricity = sqrt(1 - (b/a)^2);
range=a+b;     
% sum axes
if isnan(range); range=99;end

%% Display the eccentricity
if ploton
    fprintf('The eccentricity of the path is: %f\n', eccentricity);
    
    % Plot the fitted ellipse
    theta = linspace(0, 2*pi, 100);
    ellipse_x = a * cos(theta);
    ellipse_y = b * sin(theta);
    
    % Rotate the ellipse back to the original coordinate system
    ellipse_points = V * [ellipse_x; ellipse_y];
    
    % Translate the ellipse back to the original center
    ellipse_points(1, :) = ellipse_points(1, :) + mean_data(1);
    ellipse_points(2, :) = ellipse_points(2, :) + mean_data(2);
    
    % Plot the ellipse
    hold on;
    plot(ellipse_points(1, :), ellipse_points(2, :), 'r--');
    legend('Path', 'Fitted Ellipse');
    t=strcat('eccentricity= ', num2str(eccentricity),'     range =',num2str(range));
    title(t);
    if ~isnan(eccentricity)
      pause
      close all
    end
  end % plot on
return