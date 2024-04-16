function [MatrixOutput] = BreuerMap(MatrixInput,R,phi)
% Hard coded for 8 x 8
% To do : extend to any size


Atemp = MatrixInput([1,2,3,4],[1,2,3,4]) ; 
Btemp = MatrixInput([1,2,3,4],[5,6,7,8]) ;
Ctemp = MatrixInput([5,6,7,8],[1,2,3,4]) ; 
Dtemp = MatrixInput([5,6,7,8],[5,6,7,8]) ; 


A = subBreuer(Atemp,R,phi) ;
B = subBreuer(Btemp,R,phi) ;
C = subBreuer(Ctemp,R,phi) ;
D = subBreuer(Dtemp,R,phi) ;

MatrixOutput = [A,B ; C, D] ;

end