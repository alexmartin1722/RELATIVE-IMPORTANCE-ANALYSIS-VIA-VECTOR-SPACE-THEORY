% Author: Alex Martin
%
% Date: 7/11/22



function [Q,R] = householderQR(X)
%householderQR does the householder transformation.
%Compute the QR factorization of m,n matrix X
% Inputs:
%   X: an m,n matrix 
% Output:
%   Q: Orthogonal matrix 
%   R: Transformed matrix 

[m,n] = size(X);
Q = eye(m); %Initial orthogonal matrix 
R = X; 

for i = 1:n
    normr = norm(R(i:end,i));
    s = -sign(R(i,i));
    u = R(i,i) - s*normr;
    v = R(i:end,i)/u;
    v(1) = 1;
    tau = -s*u/normr;


    R(i:end,:) = R(i:end,:)-(tau*v)*(v'*R(i:end,:));
    Q(:,i:end) = Q(:,i:end)-(Q(:,i:end)*v)*(tau*v)';
end
Q
R

    





