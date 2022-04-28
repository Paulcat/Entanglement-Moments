function M= sigma(alpha,phi,eta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
alpha=-alpha;
M=[exp(-eta*abs(alpha*exp(1i*phi))^2),-eta*conj(alpha*exp(1i*phi))*exp(-eta*abs(alpha*exp(1i*phi))^2);
   -eta*alpha*exp(1i*phi)*exp(-eta*abs(alpha*exp(1i*phi))^2), (1-eta+((eta)^2)*abs(alpha*exp(1i*phi))^2)*exp(-eta*abs(alpha*exp(1i*phi))^2)];
M=(2*M-eye(2,2));
end
