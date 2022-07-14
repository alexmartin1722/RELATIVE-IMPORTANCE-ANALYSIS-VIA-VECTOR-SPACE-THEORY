% Author: Alex Martin
%
% Date: 7/14/22 (finished passes Q'*Q=I)
%https://www.math.usm.edu/lambers/mat610/sum10/lecture9.pdf
function [Q,R] = givens_rotationQR(X)
%givens_rotation does the Gives Rotation process.
%Compute the QR factorization of m,n matrix X
% Inputs:
%   X: an m,n matrix 
% Output:
%   Q: Orthogonal matrix 
%   R: Upper triangular  matrix 
% X = readmatrix('women.xlsx');
% X(isnan(X))=0;
[m,n] = size(X);
Q = eye(m);
R = X; 
for i = 1:n
    for j = m:-1:i+1
        G = eye(m);
        [c,s] = givens_rotation(R(j-1,i),R(j,i));
        G([j-1, j],[j-1, j]) = [c -s; s c];
        R = G'*R;
        Q = Q*G;
    end
end
R = triu(R);
R = R(1:n,1:n); %trim R to an n,n matrix 
Q = Q(1:m,1:n); %trim Q to match size of input matrix
end 

function [c,s] = givens_rotation(r1, r2)
if r2 == 0
    c = 1;
    s = 0;
else
    if abs(r2) > abs(r1)
        r = r1 / r2;
        s = 1 / sqrt(1 + r^2);
        c = s*r;
    else
        r = r2 / r1;
        c = 1 / sqrt(1 + r^2);
        s = c*r;
    end
end
end


