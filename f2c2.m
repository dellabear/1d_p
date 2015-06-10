function uc = f2c2( u )
%
% function uc = f2c2( u )
%
% Transfers a fine grid to a coarse grid by full weighting
%
% Input   
%   u     fine-grid function   
%
% Returns
%   uc    coarse-grid function

[nx] = length( u );
nxc = (nx-1)/2;
uc = zeros(nxc,1);

% do a full weighting to get value for each point on the coarse grid

  for i=1:nxc
    uc(i) = 0.5*u(2*i)+0.25*(u(2*i-1)+u(2*i+1));
  end
  
end


