% PHYSICAL DATA OF THE PROBLEM
clear problem_data

nrb = nrb4surf ([-0.5, -0.5], [0.5, -0.5], [-0.5, 0.5], [0.5 0.5]);
nrb = nrbdegelev (nrb, [1 1]);
nrb.coefs (3,:,:) = [0 0.5 0; -0.5 0 -0.5; 0 0.5 0];
problem_data.geo_name = nrb;

problem_data.drchlt_sides = 1;
problem_data.drchlt_components = {[1 2 3]};
problem_data.rotation_sides = 1;

% Physical parameters
E = 2e11;
nu = 0.3;
thickness = 0.01;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) -8000*thickness*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

% DISCRETIZATION PARAMETERS
clear method_data
deg = 3;
method_data.degree     = deg*[1 1];      % Degree of the splines
method_data.regularity = (deg-1)*[1 1];  % Regularity of the splines
method_data.nsub       = 16*[1 1];       % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad      = (deg+1)*[1 1];  % Points for the Gaussian quadrature rule
method_data.penalty_coeff = 10;

[geometry, msh, space, u] = solve_kirchhoff_love_shell (problem_data, method_data);
pts = {[1], [0.5]};
displ = sp_eval (u, space, geometry, pts)

output_file = 'hyperboloid_1patch';
vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'Displacement', 'value')
