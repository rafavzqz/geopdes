% EX_KL_SHELL_SCORDELIS_LO_ROOF: solve the Kirchhoff-Lovel shell model for the Scordelis Lo roof.

% Physical domain, defined as NURBS map
clear problem_data
radius = 25.0;
theta = 2.0 / 9.0 * pi;

nrb = nrbreverse (nrbcirc(25, [0 0 0], pi/2 - theta, pi/2 + theta));
srf = nrbextrude (nrbtform(nrb,vecrotx(pi/2)), [0 50 0]);
problem_data.geo_name = srf;

% Type of boundary conditions for each side of the domain
% Only homogeneous Dirichlet conditions have been implemented so far.
problem_data.drchlt_sides = [3 4];
problem_data.drchlt_components = {[1 3] [1 3]};

% Physical parameters
E = 4.32e8;
nu = 0.0;
thickness = 0.25;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) -90*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

% Discretization parameters
deg = 3;

clear method_data
method_data.degree     = [deg deg];
method_data.regularity = [deg-1 deg-1];
method_data.nsub       = [16 16];
method_data.nquad      = [deg+1 deg+1];

% Call to solver
[geometry, msh, space, u] = solve_kirchhoff_love_shell (problem_data, method_data);

% Postprocessing
output_file = strcat('ScordelisLo_deg',num2str(deg));
vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'Displacement')

pts{1} = 0.0; pts{2} = 0.5;
[eu,~] = sp_eval(u,space, geometry,pts);
reference_displacement = -0.3005925;
fprintf('Numerical displacement at edge: %d. \n', eu(3));
fprintf('Reference displacement at edge: %d. \n', reference_displacement);
