% Author: Alex Martin
%
% Date: 7/13/22



function [Q,R] = householderQR2(X)
%householderQR does the householder transformation.
%Compute the QR factorization of m,n matrix X
% Inputs:
%   X: an m,n matrix 
% Output:
%   Q: Orthogonal matrix 
%   R: Upper triangular  matrix 

[m,n] = size(X);
Q = eye(m);
R = X;

for i = 1:m-1
    x = R(i:m, i);
    v = -sign(x(1))*norm(x)*eye(m-i+1,1)-x;
    if norm(v)>0
        v=v/norm(v);
        P=eye(m); 
        P(i:m,i:m)=P(i:m,i:m)-2*v*v';
        R=P*R;
        Q=Q*P;
    end
end
Q  = Q * -1
R = triu(R) * -1