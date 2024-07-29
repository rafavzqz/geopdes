% MP_SOLVE_LAPLACE_C1: solve the Laplacian problem in a multipatch geometry, with C1 continuity.
%
% Example to solve the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% where \Omega is an analysis-suitable G1 multipatch domain.
%
% USAGE:
%
%  [geometry, msh, space, u] = 
%          mp_solve_laplace_C1 (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%
% OUTPUT:
%
%  geometry: array of geometry structures (see geo_load)
%  msh:      multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  space:    multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch_C1)
%  u:        the computed degrees of freedom
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011, 2013, 2015, 2017, 2022 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, space, u] = ...
              mp_solve_laplace_C1 (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch
% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
end

[edges, vertices] = vertices_struct (geometry, interfaces, boundaries, boundary_interfaces);
msh = msh_multipatch (msh, boundaries);
space = sp_multipatch_C1 (sp, msh, geometry, edges, vertices);
clear sp

% Compute and assemble the matrices 
stiff_mat = op_gradu_gradv_mp (space, space, msh, c_diff);
rhs = op_f_v_mp (space, msh, f);

% Apply Neumann boundary conditions
for iref = nmnn_sides
  gref = @(varargin) g(varargin{:}, iref);
  for bnd_side = 1:msh.boundaries(iref).nsides
    iptc = msh.boundaries(iref).patches(bnd_side);
    iside = msh.boundaries(iref).faces(bnd_side);

    msh_side = msh.msh_patch{iptc}.boundary(iside);
    sp_side = space.sp_patch{iptc}.boundary(iside);
    rhs_nmnn = op_f_v_tp (sp_side, msh_side, gref);
    [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
    rhs(Cpatch_cols) = rhs(Cpatch_cols) + Cpatch(sp_side.dofs,:).' * rhs_nmnn;
%     rhs(space.Cpatch_cols{iptc}) = rhs(space.Cpatch_cols{iptc}) + ...
%       space.Cpatch{iptc}(sp_side.dofs,:).' * rhs_nmnn;
  end
end

% Apply Dirichlet boundary conditions in weak form, by Nitsche's method
if (exist ('weak_drchlt_sides', 'var'))
  [N_mat, N_rhs] = sp_weak_drchlt_bc_laplace (space, msh, weak_drchlt_sides, h, c_diff, Cpen);
  stiff_mat = stiff_mat - N_mat;
  rhs = rhs + N_rhs;
end

% Solve the linear system
u = stiff_mat \ rhs;

end
