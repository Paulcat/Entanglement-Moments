function [MatrixOutput] = vMat(MatrixInput,R,phi)
%

% Orthogonal matrix
% We chose a specific one ; permutation of coordinate axes
% To do : check if other give different results


D = zeros(4) ; 
D(1,2) = exp(i*phi(1)) ;
D(2,1) = -exp(i*phi(1)) ;

D(3,4) = exp(i*phi(2)) ;
D(4,3) = -exp(i*phi(2)) ;

U = R*D*transpose(R) ; 

MatrixOutput = U*transpose(MatrixInput)*ctranspose(U) ;

end