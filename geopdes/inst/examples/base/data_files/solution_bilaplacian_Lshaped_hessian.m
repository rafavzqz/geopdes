function hessian = solution_bilaplacian_Lshaped_hessian (x, y)
hessian = zeros(2,2,size(x,1), size(x,2));
[th, r] = cart2pol (x, y);
th = (th < 0).*(2*acos(-1) + th) + (th >= 0) .* th;
z = 0.544483736782464;

C1 = 1/(z-1)*sin(3*pi/2*(z-1)) - 1/(z+1)*sin(3*pi/2*(z+1));
C2 = cos(3*pi/2*(z-1)) - cos(3*pi/2*(z+1));
F1 = cos((z-1)*th) - cos((z+1)*th);
F2 = 1/(z-1)*sin((z-1)*th) - 1/(z+1)*sin((z+1)*th);

psi = C1*F1 - C2*F2;  %C1 and C2 are constants
Dpsi = C1*(-(z-1)*sin((z-1)*th)+(z+1)*sin((z+1)*th)) - C2*(cos((z-1)*th)-cos((z+1)*th));
D2psi = C1*(-((z-1)^2)*cos((z-1)*th)+((z+1)^2)*cos((z+1)*th)) - C2*(-(z-1)*sin((z-1)*th)+(z+1)*sin((z+1)*th));

ur=(z+1)*r.^(z).*psi;
uth=r.^(z+1).*Dpsi;
urr=z*(z+1)*r.^(z-1).*psi;
uthr=(z+1)*r.^(z).*Dpsi;
uthth=r.^(z+1).*D2psi;

uxr=cos(th).*urr-sin(th).*(r.*uthr-uth)./(r.^2);
uxth=-sin(th).*ur+cos(th).*uthr-cos(th).*uth./r-sin(th).*uthth./r;
uyr=sin(th).*urr+cos(th).*(r.*uthr-uth)./(r.^2);
uyth=cos(th).*ur+sin(th).*uthr-sin(th).*uth./r+cos(th).*uthth./r;

uxx=cos(th).*uxr-sin(th).*uxth./r;
uxy=sin(th).*uxr+cos(th).*uxth./r;
uyy=sin(th).*uyr+cos(th).*uyth./r;

hessian(1,1,:,:) = uxx; 
hessian(1,2,:,:) = uxy; 
hessian(2,1,:,:) = uxy; 
hessian(2,2,:,:) = uyy;

end