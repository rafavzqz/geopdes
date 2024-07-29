% TEST_PLATE_MIXED_BC_G_NMNN: data function for Neumann boundary condition.

function g = bilaplacian_Lshaped_g_nmnn_8patches (x, y, ind)
  [th, r] = cart2pol (x, y);
  th = (th < 0).*(2*acos(-1) + th) + (th >= 0) .* th;
% f = @(z) sin(3*pi/2*z).^2 - z.^2 * sin(3*pi/2)^2;
  z = 0.544483736782464; %z = fsolve(f, 0.5) The one on the left is more precise, computed with bisection method

  C1 = 1/(z-1)*sin(3*pi/2*(z-1)) - 1/(z+1)*sin(3*pi/2*(z+1));
  C2 = cos(3*pi/2*(z-1)) - cos(3*pi/2*(z+1));
  F1 = cos((z-1)*th) - cos((z+1)*th);
  F2 = 1/(z-1)*sin((z-1)*th) - 1/(z+1)*sin((z+1)*th);
  DF1 = -(z-1)*sin((z-1)*th)+(z+1)*sin((z+1)*th);
  DF2 = cos((z-1)*th)-cos((z+1)*th);

  psi = C1*F1 - C2*F2;  %C1 and C2 are constants
  Dpsi = C1*DF1 - C2*DF2;

  switch (ind)
    case {4} %n=(-1,0)
      g = -((z+1)*r.^(z).* psi.*cos(th) - r.^z.*Dpsi.*sin(th));
    case {2,6} %n=(1,0)
      g = (z+1)*r.^(z).* psi.*cos(th) - r.^z.*Dpsi.*sin(th);
    case {1,3} %n=(0,-1)
      g = -((z+1)*r.^(z).* psi.*sin(th) + r.^z.*Dpsi.*cos(th));
    case {5} %n=(0,1)
      g = (z+1)*r.^(z).* psi.*sin(th) + r.^z.*Dpsi.*cos(th);
        
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end

