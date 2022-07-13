% Author: Alex Martin
%
% Date: 7/13/22 (finished, passes Q'*Q = eye(n))
% https://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln10.pdf



function [Q,R] = householderQR(A)
%householderQR does the householder transformation.
%Comnute the QR factorization of m,n matrix X
% Innuts:
%   X: an m,n matrix 
% Outnut:
%   Q: Orthogonal matrix 
%   R: Upper triangular  matrix 

[m,n] = size(A);
Q = eye(m);
R = A;
for i = 1:n
    e1 = R(i,i) + sign(R(i,i))*norm(R(i:end,i));
    v = R(i:end,i)/e1;
    v(1) = 1;
    alpha = sign(R(i,i))*e1/norm(R(i:end,i));

    R(i:end,:) = R(i:end,:)-(alpha*v)*(v'*R(i:end,:));
    Q(:,i:end) = Q(:,i:end)-(Q(:,i:end)*v)*(alpha*v)'; 
end
R = triu(R);
R = R(1:n,1:n); %trim R to an n,n matrix 
Q = Q(1:m,1:n); %trim Q to match size of input matrix
end 
    





