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
%   R: Upper triangular  matrix 

[m,n] = size(X);
Q = eye(m);
R = X;

for i = 1:n
    z = R(i,i)-(-sign(R(i,i))*norm(R(i:end,i)));
    v = R(i:end,i)/z;
    v(1) = 1;
    p = (sign(R(i,i))*z)/norm(R(i:end,i));

    R(i:end,:) = R(i:end,:)-(p*v)*(v'*R(i:end,:));
    Q(:,i:end) = Q(:,i:end)-(Q(:,i:end)*v)*(p*v)';
end
Q = Q * -1;
R = triu(R) * -1;


    





