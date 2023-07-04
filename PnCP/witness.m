% Define the normal vector
V = [a, b, c];    % Normal vector

% Create a meshgrid of points within a certain range to form a grid for the plane
x=linspace(0, 2.25, 100);
y=linspace(0, 2.25, 100);
[X, Y] = meshgrid(x, y);

% Calculate the z-coordinate of each point on the grid using the equation of the plane
Z = (d-V(1)*X - V(2)*Y)/V(3);


% Plot the plane using the meshgrid
figure;
surf(X, Y, Z);
 
 
% Set the range of the axes
axis([0, 2.25, 0, 2.25, 0, 2.25]);

hold;

plot(h6, 'FaceColor', 'red', 'MarkerFaceColor','red', 'FaceAlpha', 0.5, 'LineWidth',0.8); %h6 is the pncp polytope.
axis([0, 2.25, 0, 2.25, 0, 2.25]);
 
%scatter3(1.5,1,1, 70, '.', 'MarkerEdgeColor','black'); %the three coordinates are the point we wanted to detect.
%axis([0, 2.25, 0, 2.25, 0, 2.25]);

