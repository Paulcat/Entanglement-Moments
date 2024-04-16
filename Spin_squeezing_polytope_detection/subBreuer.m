function [MatrixOutput] = subBreuer(MatrixInput,R, phi)
% Applying the Breuer to the 4x4 subpart
sza = size(MatrixInput) ;

if isequal(sza(1),sza(2))
   if isequal(sza(1),4)
   else
       disp('Size is not 4x4 !')
   end
else
   disp('Not square matrix')
end

A = MatrixInput ;

MatrixOutput = eye(4)*trace(A) - A - vMat(A,R,phi) ;

end