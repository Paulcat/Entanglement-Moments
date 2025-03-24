clear all;

d = 3; % Matrix dimension
num_terms = 4; % Number of Hermitian terms in the sum

% Create a product manifold of Hermitian matrices (as a structure array)
manifold_list = struct(); % Initialize as a structure
for k = 1:num_terms * 2
    field_name = ['H', num2str(k)]; % Create a field name like 'H1', 'H2', etc.
    manifold_list.(field_name) = sympositivedefinitefactory(d); % Hermitian matrix manifold
end
manifold = productmanifold(manifold_list); % Now it's a structure array

% Define the cost function (linked to the SDP)
problem.M = manifold;
problem.cost = @(H) -Distinguishing_PT_PhiPos_Sets_From_Unitaries(build_W_from_H(H), d);

% Run the optimization
H0 = manifold.rand(); % Random initialization
[H_opt, fval, info] = conjugategradient(problem, H0);

% Display the optimal value
disp("Optimal value reached:");
disp(-fval);
