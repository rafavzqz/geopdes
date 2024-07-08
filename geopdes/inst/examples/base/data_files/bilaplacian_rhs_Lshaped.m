function uex = bilaplacian_rhs_Lshaped (x, y)
[th, r] = cart2pol (x, y);
th = (th < 0).*(2*acos(-1) + th) + (th >= 0) .* th;
z = 0.544483736782464;

C1 = 1/(z-1)*sin(3*pi/2*(z-1)) - 1/(z+1)*sin(3*pi/2*(z+1));
C2 = cos(3*pi/2*(z-1)) - cos(3*pi/2*(z+1));
F1 = cos((z-1)*th) - cos((z+1)*th);
F1_2 = -(z-1)^2*cos((z-1)*th) + (z+1)^2*cos((z+1)*th); 
F1_4 = (z-1)^4*cos((z-1)*th) - (z+1)^4*cos((z+1)*th); 
F2 = 1/(z-1)*sin((z-1)*th) - 1/(z+1)*sin((z+1)*th);
F2_2 = -(z-1)*sin((z-1)*th) + (z+1)*sin((z+1)*th); 
F2_4 = (z-1)^3*sin((z-1)*th) - (z+1)^3*sin((z+1)*th); 

psi = C1*F1 - C2*F2;
psi2 = C1*F1_2 - C2*F2_2;
psi4 = C1*F1_4 - C2*F2_4;

uex = (r.^(z-3)).*(((z-1)^2)*(((z+1)^2)*psi+psi2)+((z+1)^2)*psi2+psi4);

end