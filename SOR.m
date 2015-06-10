function [ u] = SOR(A,b,u)
% SOR smoother for 1D possion 
%   A: 
%   b:
%   u: initial function

n = length(b);
I = speye(n);
D = diag(diag(A));
L = tril(A)-D;

M = I-(D+1.8*L)\A;
c = (D+1.8*L)\b;
u = M*u+c;


end

