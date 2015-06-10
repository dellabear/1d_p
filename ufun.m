function [ y ] = ufun( n)
%  exact solution to the PDE
%  x: input


h = 1/(n+1);
x = (h:h:(1-h))';

for i=1:n
    
if x(i)>=0 && x(i)<=0.5
    y(i)=-x(i);
else
    y(i)=-1+x(i);
end
end

end

