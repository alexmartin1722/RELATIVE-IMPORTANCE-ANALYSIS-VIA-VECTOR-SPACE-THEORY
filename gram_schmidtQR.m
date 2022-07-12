% Author: Alex Martin
%
% Date: 7/13/22



function [Q,R] = gram_schmidtQR(X)
%gram_schmidtQR does the classic gram-schmidt process.
%Compute the QR factorization of m,n matrix X
% Inputs:
%   X: an m,n matrix 
% Output:
%   Q: Orthogonal matrix 
%   R: Upper triangular  matrix 
% http://web.mit.edu/18.06/www/Essays/gramschmidtmat.pdf
[m,n] = size(X);
Q = zeros(m,n);
R = zeros(m,n);

for i = 1:n
    v = X(:,i);
    for j = 1:i-1
        R(j,i) = Q(:j)'*X(:,i);
        v = v-R(j,i)*Q(:,j);
    end
    R(i,i) = norm(v);
    Q(:,i) = v/(R(i,i));
end