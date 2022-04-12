function [a,adag] = quantop(n)
%QUANTOP creation and annihilation operators

adag	= circshift(diag(sqrt(0:n-1)),-1);
a		= adag';

end