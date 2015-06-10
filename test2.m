%The test for multigrid solver for 1D poission equation


levels = 8;                                          % size of problem
nu1 = 5;                           % number of presmoothing iterations
nu2 = 5;                          % number of postsmoothing iterations
gamma = 1;  % number of coarse grid iterations (1=V-cycle, 2=W-cycle)
%---------------------------------------------------------------------
n = 2^(levels+2)-1;                            % number of grid points
h = 1/(n+1);
x = (h:h:(1-h))';

% function f(x)
epsilon=h;
f=zeros(size(x));

for i=1:n
if abs((x(i)-0.5)/epsilon)<1
f(i)= 2*(abs((x(i)-0.5)/epsilon)+1)/epsilon;
else
f(i)=0;
end
end

A = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);
b = f*h^2;
uc = A\b;

