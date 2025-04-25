function [distance] = Finding_coefficicents_from_Gellmann_Basis(alpha,rho)

d=3;
GellMann = zeros(d, d, d*d);

GellMann(:,:,1) = [0, 1, 0; 1, 0, 0; 0, 0, 0];                % lambda1
GellMann(:,:,2) = [0, -1i, 0; 1i, 0, 0; 0, 0, 0];             % lambda2
GellMann(:,:,3) = [1, 0, 0; 0, -1, 0; 0, 0, 0];               % lambda3
GellMann(:,:,4) = [0, 0, 1; 0, 0, 0; 1, 0, 0];                % lambda4
GellMann(:,:,5) = [0, 0, -1i; 0, 0, 0; 1i, 0, 0];             % lambda5
GellMann(:,:,6) = [0, 0, 0; 0, 0, 1; 0, 1, 0];                % lambda6
GellMann(:,:,7) = [0, 0, 0; 0, 0, -1i; 0, 1i, 0];             % lambda7
GellMann(:,:,8) = (1/sqrt(3)) * [1, 0, 0; 0, 1, 0; 0, 0, -2]; % lambda8
GellMann(:,:,9) = eye(3);                                     % 1|

% Build the linear combination
% Unitary=zeros(3);
% for i = 1:d*d
%     Unitary = Unitary + alpha(i) * GellMann(:,:,i);
% end
% 
% Witness=kron(Unitary,Unitary);
full_vector = zeros(81, 1); % Step 1: Create a zero vector (81x1)
indices = randperm(81, 18); % Step 3: Select 18 unique random positions

full_vector(indices) = alpha; % Step 4: Insert values at random positions
Witness = zeros(d*d,d*d);

for i = 1:d*d
    for j=1:d*d
    Witness = Witness + full_vector(i*j)*kron(GellMann(:,:,i),GellMann(:,:,j));
    end
end



distance = norm(Witness - rho);

end

