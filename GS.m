function [ u] = GS(A,b,u)
% Gauss-Seidel smoother for 1D possion 
%   A: 
%   b:
%   u: initial function

n = length(b);
I = speye(n);
D = diag(diag(A));
L = tril(A)-D;

M = I-(D+L)\A;
c = (D+L)\b;
u = M*u+c;


end

