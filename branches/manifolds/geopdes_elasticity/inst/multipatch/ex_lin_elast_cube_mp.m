% EX_LIN_ELAST_CUBE_MP: solve linear elasticity problem in a multipatch cube.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_2cubesa.txt';
% You can see how the patches should match in the files
%  geo_2cubesa{b,c,d,e,f,g,h}.txt
% In all the eight cases the result must be the same

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
E  =  1; nu = 0.3; 
problem_data.lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x))); 
problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
fx = @(x, y, z) problem_data.mu_lame(x, y, z) .* (2*cos(x) - 4*x.*y.^2.*z) + ...
                problem_data.lambda_lame(x, y, z) .* (cos(x) - 4*x.*y.^2.*z);
fy = @(x, y, z) problem_data.mu_lame(x, y, z) .* (2*sin(y).*z - 4*x.^2.*y.*z) + ...
                problem_data.lambda_lame(x, y, z) .* (z.*sin(y) - 4*x.^2.*y.*z);
fz = @(x, y, z) problem_data.mu_lame(x, y, z) .* (-2*(y.*z).^2 - 2*(x.*z).^2 - cos(y) - 4*(x.*y).^2) - ...
                problem_data.lambda_lame(x, y, z) .* (2*(x.*y).^2 + cos(y));
problem_data.f = @(x, y, z) cat(1, ...
                    reshape (fx (x,y,z), [1, size(x)]), ...
                    reshape (fy (x,y,z), [1, size(x)]), ...
                    reshape (fz (x,y,z), [1, size(x)]));
problem_data.h = @(x, y, z, ind) cat(1, ...
                    reshape (cos(x), [1, size(x)]), ...
                    reshape (sin(y).*z, [1, size(x)]), ...
                    reshape ((x.*y.*z).^2, [1, size(x)]));

% Exact solution (optional)
uxex = @(x, y, z) cos(x);
uyex = @(x, y, z) sin(y).*z;
uzex = @(x, y, z) (x.*y.*z).^2;
problem_data.uex  = @(x, y, z) cat(1, ...
		      reshape (uxex (x,y,z), [1, size(x)]), ...
		      reshape (uyex (x,y,z), [1, size(x)]), ...
		      reshape (uzex (x,y,z), [1, size(x)]));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [2 2 2];     % Degree of the bsplines
method_data.regularity = [1 1 1];     % Regularity of the splines
method_data.nsub       = [2 2 2];     % Number of subdivisions
method_data.nquad      = [3 3 3];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u, gnum] = ...
              mp_solve_linear_elasticity_3d (problem_data, method_data);

% 4) POST-PROCESSING. 
% 4.1) EXPORT TO PARAVIEW
output_file = 'lin_elast_cube_mp_Deg2_Reg1_Sub2';

vtk_pts = {linspace(0, 1, 10), linspace(0, 1, 10), linspace(0, 1, 10)};
fprintf ('results being saved in: %s_displacement.pvd\n \n', output_file)
mp_sp_to_vtk (u, space, geometry, gnum, vtk_pts, sprintf ('%s_displacement', output_file), 'displacement')

% 4.2) COMPARISON WITH THE EXACT SOLUTION
npatch = numel (geometry);
for iptc = 1:npatch
  error_l2(iptc) = ...
        sp_l2_error (space{iptc}, msh{iptc}, u(gnum{iptc}), problem_data.uex);
end
error_l2 = sqrt (sum (error_l2 .* error_l2))

%!demo
%! ex_lin_elast_cube_mp
