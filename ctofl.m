function u = ctofl( uc )
%
% function u = ctof( uc )
%
% Transfers a coarse grid to a fine grid
%
% Input
%   uc    coarse-grid function
%
% Returns
%   u     fine-grid function

[nx] = length( uc );
nxf = 2*nx+1;
u = zeros(nxf,1);

% transfer by copying where the grids line up

  for i=1:nx
    u(2*i) = uc(i);
  end

%quaratic interpolation in x

  for i=1:2:nxf
      if i==1
         u(i)=(-1/8)*u(i+3)+3/4*u(i+1)+0;         
      elseif i==nxf         
         u(i)=(-1/8)*u(i-3)+3/4*u(i-1)+0;
      elseif i==nxf-2
         u(i)=(-1/8)*u(i-3)+3/4*u(i-1)+3/8*u(i+1);
      else
         u(i)=(-1/8)*u(i+3)+3/4*u(i+1)+3/8*u(i-1);
      end
  end