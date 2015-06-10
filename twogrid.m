function [x1,x2] = twogrid(A,b,nu1,nu2,gamma,x0)
%TWOGRID
%    Recursive twogrid cycle for 1d Poisson problem.
%    nu1 = number of presmoothing iterations (Gauss-Seidel)
%    nu2 = number of postsmoothing iterations (Gauss-Seidel)
%    gamma = number of coarse grid iterations (1=V-cycle, 2=W-cycle)
%    x0 = starting vector (0 if not prescribed)

global level t

n = length(b);
fs=n+1;
levels=log(n+1)/log(2)-3;
h = 1/(n+1);
xx = (h:h:(1-h))';

x_c = A\b;

if level<1
   x = A\b;
   x1=x;x2=x;
   % solve exactly at coarsest level
else  
   I = spdiags(ones(n-2,1)*[1 2 1],-2:0,n,n-2); % create interpolation
   I = I(:,1:2:end)/2; R = I'/2;            % and restriction matrices
   if nargin<6, x = b*0; 
   else x = x0; end           % starting vector
   
   for i = 1:nu1, x = GS(A,b,x); end                       % presmoothing
   
   E1=abs(x-x_c);
   subplot(2,2,1) 
   plot(xx, E1,'co','Linewidth',2);
   hold on
   str = sprintf('Error value after smoothing at level %d', levels);
   title(str);
   hold off
   
        ws = 2*pi/n;
        wnorm = -pi:ws:pi;
    	wnorm = wnorm(1:length(x));
        w = wnorm*fs;
   
   F1=abs(fftshift(fft(E1)));
   subplot(2,2,2) 
   plot(w,F1,'co','Linewidth',2);
   str = sprintf('Error Frequency after smoothing at level %d', levels);
   title(str);
   hold off
   
   r = b-A*x;                                       % compute residual
   
   % f2c1 is ENO based weighting
   % f2c2 is normal full weighting
   %rh = f2c1(r);                        % restrict residual to coarse grid
   %rh = f2c2(r);
      rh = R*r;  
   
   level = level-1; 

   %plot([t-1 t],[level+1 level],'bo-')
   
   eh = rh*0;                                        % starting vector
   
   for i = 1:gamma
      [eh1,eh2] = twogrid(R*A*I,rh,nu1,nu2,gamma,eh); % coarse grid iteration
   end
   
 
   
   % interpolate error
      t = t + 1;    
   e1 = ctofo(eh1);    %e1 is for ENO interpolation
   e2 = ctofl(eh2);    %e2 is for normal quaratic interpolation

%    if t <=2
%        e1=e2;
%    end
%    
   
   %t = t+1; level = level+1; 
   %plot([t-1 t],[level-1 level],'bo-')
   
   % update solution
   x1 = x + e1;
   x2 = x + e2;  

   E0 = x-x_c;
   
   for i = 1:nu2, x1 = GS(A,b,x1);x2 = GS(A,b,x2); 
   E1=abs(x1-x_c);
   E2=abs(x2-x_c);
   subplot(2,2,3) 
   plot(xx,E1,'co','Linewidth',2);
   hold on
   plot(xx,E2,'ro','Linewidth',2);
   str = sprintf('Error value during postsmoothing at level %d', levels);
   title(str);
   legend('ENO','Tradiational');
   hold off
   
   
        ws = 2*pi/n;
        wnorm = -pi:ws:pi;
    	wnorm = wnorm(1:length(x));
        w = wnorm*fs;
   
   
   F1=abs(fftshift(fft(E1)));
   F2=abs(fftshift(fft(E2)));
   subplot(2,2,4) 
   plot(w,F1,'co','Linewidth',2);
   hold on
   plot(w,F2,'ro','Linewidth',2); 
   legend('ENO','Tradiational');
     
   str = sprintf('Error frequency during postsmoothing at level %d', levels);
   title(str);
   hold off
  
   end                      % postsmoothing
%    if n<=20
%     x1=x2;
%    end
end