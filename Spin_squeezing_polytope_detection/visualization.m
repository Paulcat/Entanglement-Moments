% A script to run them all
% To DO : put it into functions

% Parameters
epsilon = 0.01 ;
N = 4;

% ALgebraic
h1 = algebraic(N) ; 

disp('Plot of algebraic polytope from Spin Squeezing inequalities')

plot(h1, 'FaceColor', 'cyan', 'MarkerFaceColor','blue', 'FaceAlpha', 0.5, 'LineWidth',0.8);
hold ; 

% Physical
[PR, NPR] = PhysRegPoints1(N) ; 
disp('Computation of physical points finished')

h2 = PhysRegPointsBoxes2(N,PR,NPR); 
disp('Improvement of the computation of physical points completed')
% physical points
disp('Plot of physical points')
plot(h2, 'FaceColor', 'white', 'MarkerFaceColor','grey', 'FaceAlpha', 0.5, 'LineWidth',0.8);

% ppt physical points
[h3,h4] = pptNqubits3(N) ; 
% plot(h3, 'FaceColor', 'red', 'MarkerFaceColor','yellow', 'FaceAlpha', 0.5, 'LineWidth',0.8);
% 
% %  NPT physical points
% plot(h4, 'FaceColor', 'grey', 'MarkerFaceColor','black', 'FaceAlpha', 0.5, 'LineWidth',0.8);
% % 
% % 
% % boxes5 ; 