function uc = f2c1( u )
%
% function uc = f2c1( u )
%
% Transfers a fine grid to a coarse grid by choosing the oscillatory part
%
% Input   
%   u     fine-grid function   
%
% Returns
%   uc    coarse-grid function

[nx] = length( u );
nxc = (nx-1)/2;
uc = zeros(nxc,1);


%ENO interpolation in x

  for i=1:nxc
      if i==1
         p = [abs(u(2*i-1)-2*u(2*i)+u(2*i+1)), abs(u(2*i)-2*u(2*i+1)+u(2*i+2))];
         [~,b]=max(p);
         if b==1 
             uc(i)= 0.5*u(2*i)+0.25*(u(2*i-1)+u(2*i+1));
         else
             uc(i)= 0.5*u(2*i)+u(2*i+1)-0.5*u(2*i+2);
         end
         
      elseif i==nxc
         p = [abs(u(2*i-1)-2*u(2*i)+u(2*i+1)), abs(u(2*i)-2*u(2*i-1)+u(2*i-2))];
         [~,b]=max(p);
         if b==1 
             uc(i)= 0.5*u(2*i)+0.25*(u(2*i-1)+u(2*i+1));
         else
             uc(i)= 0.5*u(2*i)+u(2*i-1)-0.5*u(2*i-2);
         end
         
      else
         p = [abs(u(2*i-1)-2*u(2*i)+u(2*i+1)), abs(u(2*i)-2*u(2*i+1)+u(2*i+2)), abs(u(2*i)-2*u(2*i-1)+u(2*i-2))];
         [~,b]=max(p);
         if b==1 
             uc(i)= 0.5*u(2*i)+0.25*(u(2*i-1)+u(2*i+1));
         elseif b==2 
             uc(i)= 0.5*u(2*i)+u(2*i+1)-0.5*u(2*i+2);
         else 
             uc(i)= 0.5*u(2*i)+u(2*i-1)-0.5*u(2*i-2);
         end
      end
  end
