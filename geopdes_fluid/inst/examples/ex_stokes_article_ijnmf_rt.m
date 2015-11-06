% EX_STOKES_ARTICLE_IJNMF_RT: solve the Stokes problem in the driven cavity with generalized Raviart-Thomas elements.
%
% This is the example of Section 4.1.2 in the paper
%      A.Buffa, C.de Falco, G. Sangalli, 
%      IsoGeometric Analysis: Stable elements for the 2D Stokes equation
%      Internat. J. Numer. Methods Fluids, 2011

% Physical domain, defined from the aspect ratio using the NURBS toolbox
aspect_ratio = 0.9;
nrb_square = nrb4surf ([0 0], [1 0], [0 1], [1 1]);
geo_name = nrbtform (nrb_square, vecscale ([1 aspect_ratio]));

% Type of boundary conditions for each side of the domain
drchlt_sides = 1:4;

% Physical parameters, and force term
viscosity = @(x, y) ones (size (x));
f  = @(x, y) zeros ([2, size(x)]);

% Discretization parameters
element_name = 'rt';   % Element type for discretization
degree       = [ 3  3];  % Degree of the splines
regularity   = [ 2  2];  % Regularity of the splines
nsub         = [10 10];  % Number of subdivisions
nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% load geometry
geometry = geo_load (geo_name);

% Compute the mesh structure using the finest mesh
[msh_breaks, der2] = msh_set_breaks (element_name, geometry.nurbs.knots, nsub);
rule               = msh_gauss_nodes (nquad);
[qn, qw]           = msh_set_quad_nodes (msh_breaks, rule);
msh                = msh_2d (msh_breaks, qn, qw, geometry, 'der2', der2);

% Compute the space structures
[space_v, space_p, PI] = sp_bspline_fluid_2d (element_name, ...
                          geometry.nurbs.knots, nsub, degree, regularity, msh);

% Assemble the matrices
A = op_gradu_gradv_tp (space_v, space_v, msh, viscosity); 
B = PI' * op_div_v_q_tp (space_v, space_p, msh);
M = op_u_v_tp (space_p, space_p, msh, @(x,y) ones (size (x))); 
E = sum (M, 1) * PI / sum (sum (M)); 
F = op_f_v_tp (space_v, msh, f);

vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

% Apply Dirichlet boundary conditions without using the L2-projection
% This is necessary to have a divergence exactly equal to zero
vel(space_v.boundary(3).comp_dofs{1}) = -1;
vel(space_v.boundary(4).comp_dofs{1}) = 1;

drchlt_dofs = [];
for iside = drchlt_sides; 
  drchlt_dofs = union (drchlt_dofs, space_v.boundary(iside).dofs);
end

int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);
rhs_dir  = -A(int_dofs, drchlt_dofs) * vel(drchlt_dofs);

% Solve the linear system
mat = [ A(int_dofs, int_dofs), -B(:,int_dofs).',               sparse(nintdofs, 1);
       -B(:,int_dofs),          sparse(size (B,1), size(B,1)), E';
       sparse(1, nintdofs),     E,                             0];
rhs = [F(int_dofs) + rhs_dir; 
       B(:, drchlt_dofs)*vel(drchlt_dofs);
       0];

sol = mat \ rhs;
vel(int_dofs) = sol(1:nintdofs);
press = PI * sol(1+nintdofs:end-1);

% Export to ParaView
output_file = 'Driven_cavity_RT_Deg3_Reg2_Sub10';

fprintf ('The result is saved in the files %s \n and %s \n \n', ...
         [output_file '_vel'], [output_file '_press']);
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
sp_to_vtk (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], 'vel')

% Plot in Matlab
[eu, F] = sp_eval (vel, space_v, geometry, vtk_pts);
[X,  Y] = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
figure ()
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title ('Computed solution')

[div, F] = sp_eval (vel, space_v, geometry, vtk_pts, 'divergence');
figure ()
surf (X, Y, div)
view (2)
axis equal
title ('Computed divergence')

%!demo
%! ex_stokes_article_ijnmf_rt
