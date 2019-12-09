% EX_KL_SHELL_SQUARE: solve the Kirchhoff-Lovel shell model in the square.

% Physical domain, defined as NURBS map
clear problem_data
base = 10;
p11 = [0 0 1]; p12 = [base 0 1]; p21 = [0 base 1]; p22 = [base base 1];
problem_data.geo_name = nrb4surf (p11,p12,p21,p22);

% Type of boundary conditions for each side of the domain
% Only homogeneous Dirichlet conditions have been implemented so far.
problem_data.drchlt_sides = [1 2 3 4];
problem_data.drchlt_components = {[2 3] [2 3] [1 3] [1 3]};

% Physical parameters
E = 1e6;
nu = 0.;
thickness = (12 * (1 - nu*nu) / E)^(1/3);

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) 4*16*pi^4/(base^4)*sin(2*pi*x/base).*sin(2*pi*y/base);

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

% Exact solution
ux = @(x, y, z) zeros(size(x));
uy = @(x, y, z) zeros(size(x));
uz = @(x, y, z) sin(2*pi*x/base).*sin(2*pi*y/base);

problem_data.uex     = @(x, y, z) cat(1, ...
    reshape (ux (x,y,z), [1, size(x)]), ...
    reshape (uy (x,y,z), [1, size(x)]), ...
    reshape (uz (x,y,z), [1, size(x)]));

problem_data.h     = @(x, y, z, iside) cat(1, ...
    reshape (zeros(size(x)), [1, size(x)]), ...
    reshape (zeros(size(x)), [1, size(x)]), ...
    reshape (zeros(size(x)), [1, size(x)]));

% components of the gradient
ux_dx = @(x, y, z) zeros(size(x));
ux_dy = @(x, y, z) zeros(size(x));
ux_dz = @(x, y, z) zeros(size(x));
uy_dx = @(x, y, z) zeros(size(x));
uy_dy = @(x, y, z) zeros(size(x));
uy_dz = @(x, y, z) zeros(size(x));
uz_dx = @(x, y, z) 2*pi/(base^1)*cos(2*pi*x/base).*sin(2*pi*y/base);
uz_dy = @(x, y, z) 2*pi/(base^1)*sin(2*pi*x/base).*cos(2*pi*y/base);
uz_dz = @(x, y, z) zeros(size(x));


problem_data.graduex = @(x, y, z) cat (1, ...
    reshape (ux_dx(x,y,z), [1, size(x)]), ...
    reshape (uy_dx(x,y,z), [1, size(x)]), ...
    reshape (uz_dx(x,y,z), [1, size(x)]), ...
    reshape (ux_dy(x,y,z), [1, size(x)]), ...
    reshape (uy_dy(x,y,z), [1, size(x)]), ...
    reshape (uz_dy(x,y,z), [1, size(x)]), ...
    reshape (ux_dz(x,y,z), [1, size(x)]), ...
    reshape (uy_dz(x,y,z), [1, size(x)]), ...
    reshape (uz_dz(x,y,z), [1, size(x)]));

% components of the hessian - u_x
ux_dxx = @(x, y, z) zeros(size(x));
ux_dxy = @(x, y, z) zeros(size(x));
ux_dxz = @(x, y, z) zeros(size(x));
ux_dyx = @(x, y, z) zeros(size(x));
ux_dyy = @(x, y, z) zeros(size(x));
ux_dyz = @(x, y, z) zeros(size(x));
ux_dzx = @(x, y, z) zeros(size(x));
ux_dzy = @(x, y, z) zeros(size(x));
ux_dzz = @(x, y, z) zeros(size(x));

% components of the hessian - u_y
uy_dxx = @(x, y, z) zeros(size(x));
uy_dxy = @(x, y, z) zeros(size(x));
uy_dxz = @(x, y, z) zeros(size(x));
uy_dyx = @(x, y, z) zeros(size(x));
uy_dyy = @(x, y, z) zeros(size(x));
uy_dyz = @(x, y, z) zeros(size(x));
uy_dzx = @(x, y, z) zeros(size(x));
uy_dzy = @(x, y, z) zeros(size(x));
uy_dzz = @(x, y, z) zeros(size(x));

% components of the hessian - u_z
uz_dxx = @(x, y, z) -4*pi^2/(base^2)*sin(2*pi*x/base).*sin(2*pi*y/base);
uz_dxy = @(x, y, z) 4*pi^2/(base^2)*cos(2*pi*x/base).*cos(2*pi*y/base);
uz_dxz = @(x, y, z) zeros(size(x));
uz_dyx = @(x, y, z) 4*pi^2/(base^2)*cos(2*pi*x/base).*cos(2*pi*y/base);
uz_dyy = @(x, y, z) -4*pi^2/(base^2)*sin(2*pi*x/base).*sin(2*pi*y/base);
uz_dyz = @(x, y, z) zeros(size(x));
uz_dzx = @(x, y, z) zeros(size(x));
uz_dzy = @(x, y, z) zeros(size(x));
uz_dzz = @(x, y, z) zeros(size(x));

problem_data.hessuex = @(x, y, z) cat (1, ...
    reshape (ux_dxx(x,y,z), [1, size(x)]), ...
    reshape (uy_dxx(x,y,z), [1, size(x)]), ...
    reshape (uz_dxx(x,y,z), [1, size(x)]), ...
    reshape (ux_dyx(x,y,z), [1, size(x)]), ...
    reshape (uy_dyx(x,y,z), [1, size(x)]), ...
    reshape (uz_dyx(x,y,z), [1, size(x)]), ...
    reshape (ux_dzx(x,y,z), [1, size(x)]), ...
    reshape (uy_dzx(x,y,z), [1, size(x)]), ...
    reshape (uz_dzx(x,y,z), [1, size(x)]), ...
    reshape (ux_dxy(x,y,z), [1, size(x)]), ...
    reshape (uy_dxy(x,y,z), [1, size(x)]), ...
    reshape (uz_dxy(x,y,z), [1, size(x)]), ...
    reshape (ux_dyy(x,y,z), [1, size(x)]), ...
    reshape (uy_dyy(x,y,z), [1, size(x)]), ...
    reshape (uz_dyy(x,y,z), [1, size(x)]), ...
    reshape (ux_dzy(x,y,z), [1, size(x)]), ...
    reshape (uy_dzy(x,y,z), [1, size(x)]), ...
    reshape (uz_dzy(x,y,z), [1, size(x)]), ...
    reshape (ux_dxz(x,y,z), [1, size(x)]), ...
    reshape (uy_dxz(x,y,z), [1, size(x)]), ...
    reshape (uz_dxz(x,y,z), [1, size(x)]), ...
    reshape (ux_dyz(x,y,z), [1, size(x)]), ...
    reshape (uy_dyz(x,y,z), [1, size(x)]), ...
    reshape (uz_dyz(x,y,z), [1, size(x)]), ...
    reshape (ux_dzz(x,y,z), [1, size(x)]), ...
    reshape (uy_dzz(x,y,z), [1, size(x)]), ...
    reshape (uz_dzz(x,y,z), [1, size(x)]));

% Discretization parameters
deg = 3;

clear method_data
method_data.degree     = [deg deg];
method_data.regularity = [deg-1 deg-1];
method_data.nsub       = [16 16];
method_data.nquad      = [deg+1 deg+1];

% Call to solver
[geometry, msh, space, u] = solve_kirchhoff_love_shell (problem_data, method_data);

% Post-processing
[~, ~, err_l2, err_h2s, err_h1s] = sp_h2_error (space, msh, u, ...
            problem_data.uex, problem_data.graduex, problem_data.hessuex);

fprintf('Error in L2 norm = %g\n', err_l2);
fprintf('Error in H1 semi norm = %g\n', err_h1s);
fprintf('Error in H2 semi norm = %g\n', err_h2s);

deformed_geometry = geo_deform (u, space, geometry);
nrbkntplot(deformed_geometry.nurbs)