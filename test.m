%The test for multigrid solver for 1D poission equation


levels = 1;                                          % size of problem
nu1 = 5;                           % number of presmoothing iterations
nu2 = 5;                          % number of postsmoothing iterations
gamma = 1;  % number of coarse grid iterations (1=V-cycle, 2=W-cycle)
%---------------------------------------------------------------------
n = 2^(levels+3)-1;                            % number of grid points
h = 1/(n+1);
x = (h:h:(1-h))';

% function f(x)
epsilon=h;
f=zeros(size(x));

% for i=1:n
% if abs((x(i)-0.5)/epsilon)<1
% f(i)= -(abs((x(i)-0.5)/epsilon)+1)/epsilon;
% % elseif abs((x(i)-0.8555)/epsilon)<1
% % f(i)= -2*(abs((x(i)-0.8555)/epsilon)+1)/epsilon;
% % elseif abs((x(i)-0.0555)/epsilon)<1
% % f(i)= -2*(abs((x(i)-0.0555)/epsilon)+1)/epsilon;
% else
% f(i)=0;
% end
% end


f = pi^2*(sin(3*pi*x)+4^2*sin(pi*4*x)+9^2*sin(pi*10*x));

A = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);
b = f*h^2;
uc = A\b;

global t level, t = 0; level = levels;
%clf, subplot(2,2,3:4), hold on



[u1,u2] = twogrid(A,b,nu1,nu2,gamma);

err1 = abs(u1-uc); %error
err2 = abs(u2-uc);
error1 = norm(err1(:),2)*h 
error2 = norm(err2(:),2)*h               %L2 norm of the error

% hold off, axis tight
% subplot(2,2,3), plot(x,u1,'c.-',x,u2,'r.-')%,x,uc,'b.--')
% title('multigrid approximation by ENO and normal interpolation')
% subplot(2,2,4), plot(x,err1,'c.-',x,err2,'r.-')
% title('Error for ENO and normal interpolation')