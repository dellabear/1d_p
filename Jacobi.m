function [ u] = Jacobi(A,b,u)
% Jacobi smoother for 1D possion 
%   A: 
%   b:
%   u: initial function

n = length(b);
I = speye(n);
D = diag(diag(A));
L = tril(A)-D;

M = I-D\A;
c = D\b;
u = M*u+c;


end

