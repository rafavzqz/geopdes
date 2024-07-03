% PHYSICAL DATA OF THE PROBLEM
clear problem_data

problem_data.geo_name = 'geo_hyperboloid_ASG1.txt';

problem_data.drchlt_sides = 2;
problem_data.drchlt_components = {[1 2 3]};
problem_data.rotation_sides = 2;

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
deg = 4;
method_data.degree     = deg*[1 1];      % Degree of the splines
method_data.regularity = (deg-2)*[1 1];  % Regularity of the splines
method_data.nsub       = [16 16];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad      = (deg+1)*[1 1];  % Points for the Gaussian quadrature rule
method_data.penalty_coeff = 10;

[geometry, msh, space, u] = ...
           mp_solve_kirchhoff_love_shell_C1_scalarspace (problem_data, method_data);
% K = op_KL_shells_mp (space, space, msh, problem_data.E_coeff, problem_data.nu_coeff, thickness);
% energy = 0.5 * u.' * K * u;
% displ = sp_eval_phys (u, space, geometry, [0.5; 0; 0.25]);
         
output_file = 'hyperboloid_3patch';
vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'Displacement', 'value')
