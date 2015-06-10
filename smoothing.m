function smoothing
clear 

%    Sets up a 1d Poisson problem with an oscillatory right hand
%    side and applies various iterative solvers to it.
%    Observe the decay rate of various error frequencies:
%    Jacobi smooths medium frequencies best, but fails for
%    grid-scale oscillations. Weighted-Jacobi remedies this
%    problem at the expense of a slower decay rate for small
%    frequencies. Gauss-Seidel performs well as a smoother.
%    SOR smooths best, but it is also more expensive in practice.
%    Observe how the error gets smoothed very quickly, but the
%    actual convergence to the correct solution is very slow.

% set up 1d Poisson problem with singular right hand side
levels = 5;     
n = 2^(levels+2)-1;
h = 1/(n+1);
x = (h:h:(1-h))';

b = 11;
sx = 2*pi*b*x.^2;
f = 1e2*sx.*sin(sx);


% function f(x)
% epsilon=h;
% f=zeros(size(x));
% 
% for i=1:n
% if abs((x(i)-0.5)/epsilon)<1
% f(i)= 2*(abs((x(i)-0.5)/epsilon)+1)/epsilon;
% else
% f(i)=0;
% end
% end

A = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);
b = f*h^2;
uc = A\b;

I = speye(n);
D = diag(diag(A));
L = tril(A)-D;
%------------------------------------------------------------------
% different preconditioners
name = {'Jacobi','weighted Jacobi','Gauss-Seidel','SOR'};
P = {D,3/2*D,D+L,D+1.8*L};
%------------------------------------------------------------------
np = length(P);
for k = 1:np
    M{k} = I-P{k}\A;
    c{k} = P{k}\b;
    u{k} = zeros(n,1);
end
for i = 1:1000
    clf
    for k = 1:np, u{k} = M{k}*u{k}+c{k}; end
    
    subplot(1,3,1)
    plot(x,uc,'k.-')
    hold all, for k = 1:np, plot(x,u{k},'.-'), end, hold off
    title(sprintf('solution after %d steps',i))
    
    subplot(1,3,2)
    hold all, 
    for k = 1:np, 
        a=abs(fft(uc-u{k},64))/64;
        a=a(1:50);
        plot(a,'.-'), 
    end, 
    hold off
    title(sprintf('error frequence %d steps',i))
    
    subplot(1,3,3)
    hold all, for k = 1:np, plot(x,uc-u{k},'.-'), end, hold off
    legend(name)
    title(sprintf('error after %d steps',i))
      drawnow
end