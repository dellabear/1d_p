function [ u] = wJacobi(A,b,u)
% Weighted-Jacobi smoother for 1D possion 
%   A: 
%   b:
%   u: initial function

n = length(b);
I = speye(n);
D = diag(diag(A));
L = tril(A)-D;

M = I-(3/2*D)\A;
c = (3/2*D)\b;
u = M*u+c;


end

