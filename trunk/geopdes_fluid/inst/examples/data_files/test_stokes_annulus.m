% TEST_STOKES_ANNULUS: data file for Stokes problem in one eighth of a ring.

degree       = [2 2];   % Degree of the spline space for the pressure
regularity   = [1 1];   % Regularity of the spline space for the pressure
nquad        = [4 4];   % Number of quadrature points
nbreaks      = [10 10]; % Number of breaks (level of refinement)

% Type of boundary conditions
drchlt_sides = 1:4;

% Geometry description
geo_name     = 'annulus.mat';

% Discrete space (function to compute it)
sp_type = 'th';
switch sp_type
  case 'th'
    test_set_th_space; 
  case 'ndl'
    test_set_ndl_space;
  case 'rt'
    test_set_rt_space;
end

% Functions to compute the viscosity and the Dirichlet boundary condition
mu = @(x, y) ones (size (x));
h            = @(x, y, iside) zeros ([2, size(x)]);

% Function to compute the right hand-side
fx = @(x, y) (4*(20*(98-57*x.^2).*y.^7+42*x.*(34*x.^2-75).*y.^6 ...
                 -3*(430*x.^4-1480*x.^2+1023).*y.^5 ...
                 +15*x.*(70*x.^4-310*x.^2+297).*y.^4 ...
                 +3*x.*(x.^4-5*x.^2+4).^2 ...
                 +2*(-290*x.^6+1500*x.^4-2079*x.^2+680).*y.^3 ...
                 +6*x.*(38*x.^6-255*x.^4+495*x.^2-260).*y.^2 ...
                 +(-75*x.^8+520*x.^6-1089*x.^4+720*x.^2-112).*y ...
                 +603*x.*y.^8-355*y.^9)-y./(x.^2+y.^2));

fy = @(x, y) (x./(x.^2+y.^2)+4*(5*x.^9-27*x.^8.*y+20*x.^7.*(15*y.^2-2) ...
                                +14*x.^6.*y.*(15-38*y.^2)+3*x.^5.*(290*y.^4-520*y.^2+33) ...
                                -15*x.^4.*y.*(70*y.^4-170*y.^2+33) ...
                                +x.^3.*(860*y.^6-3000*y.^4+2178*y.^2-80) ...
                                +18*x.^2.*y.*(-34*y.^6+155*y.^4-165*y.^2+20) ...
                                +x.*(285*y.^8-1480*y.^6+2079*y.^4-720*y.^2+16) ...
                                -67*y.^9+450*y.^7-891*y.^5+520*y.^3-48*y));
f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));


% Exact solution, to compute the errors
uxex = @(x, y) (2*(x-y).*y.*(-4+x.^2+y.^2).*(-1+x.^2+y.^2).* ...
                (x.^5-8*y-2*x.^4.*y+20*y.^3-6*y.^5+2*x.^2.*y.*(5-4*y.^2)+x.^3.*(-5 + 6*y.^2)+ ...
                 x.*(4+5*y.^2.*(-3+y.^2))));

uyex = @(x, y) (-2*(x - y).*y.^2.*(-4 + x.^2 + y.^2).*(-1 + x.^2 + y.^2).* ...
                (4 + 5*x.^2.*(-3 + x.^2) + 2*x.*(5 - 2*x.^2).*y + ...
                 (-5 + 6*x.^2).*y.^2 - 4*x.*y.^3 + y.^4));

uex = @(x, y) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));

pex = @(x, y) (-(pi/8) + atan (y./x));

graduex = @test_stokes_annulus_graduex;

% Output file for Paraview
output_file  = 'test_stokes_annulus';

% Points for post-processing
vtk_pts = {linspace(0, 1, 20)', linspace(0, 1, 20)'};

