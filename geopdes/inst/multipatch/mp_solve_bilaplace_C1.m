% MP_SOLVE_BILAPLACE_C1: solve the bilaplacian problem in a multipatch geometry.
%
% USAGE:
%
%  [geometry, msh, space, u] = 
%          mp_solve_bilaplace_C1 (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient
%    - f:            source term
%    - g:            function for Neumann condition
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
%  geometry: array of geometry structures (see mp_geo_load)
%  msh:      multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  space:    multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch_C1)
%  u:        the computed degrees of freedom
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010--2022 Rafael Vazquez
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
              mp_solve_bilaplace_C1 (problem_data, method_data)

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
% space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
space = sp_multipatch_C1 (sp, msh, geometry, edges, vertices);
clear sp

% Compute and assemble the matrices 
stiff_mat = op_laplaceu_laplacev_mp (space, space, msh, c_diff);
rhs       = op_f_v_mp (space, msh, f);

% Apply boundary conditions
u = zeros (space.ndof, 1);
if (isfield(problem_data, 'graduex') && isfield(problem_data, 'uex'))
  [u_drchlt, drchlt_dofs, add_int_dofs] = sp_bilaplacian_drchlt_C1_exact (space, msh, drchlt_sides, uex, graduex);
else
  [u_drchlt, drchlt_dofs, add_int_dofs] = sp_bilaplacian_drchlt_C1 (space, msh, drchlt_sides, h, g);
end
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
add_dofs = []; %this will contain the "boundary" vertex dofs which have been removed from drchlt_dofs
for bv = 1 : numel(add_int_dofs)
  vertex_number = add_int_dofs(bv).vertex_number;
  add_dofs = [add_dofs space.dofs_on_vertex{vertex_number}(add_int_dofs(bv).function_index)];
end
%By subtracting add_dofs (later they will be added again), int_dofs contain ony the "original" interior dofs
int_dofs = setdiff (int_dofs, add_dofs); 

%We assemble the (pieces of the) stiffness matrix, the rhs (and its correction taking 
%into account the Dirichlet conditions), and the basis change matrix (we will need it 
%to go from the basis with kernel vectors obtained when examnining the dirichlet conditions 
%to the usual basis)
add_mat = [];
add_rhs = [];
add_rb = [];
s_vec = numel(size(add_int_dofs, 1), size(add_int_dofs, 1));
B_change = speye(space.ndof); %basis change matrix

for bv = 1 : numel(add_int_dofs)
  vertex_number = add_int_dofs(bv).vertex_number;
  dofs_on_vertex = space.dofs_on_vertex{vertex_number};
  kernel_coeffs = add_int_dofs(bv).kernel_coeffs;
  B_change(dofs_on_vertex, dofs_on_vertex(add_int_dofs(bv).function_index)) = kernel_coeffs;
  
  add_mat = [add_mat; kernel_coeffs.' * stiff_mat(dofs_on_vertex, int_dofs) ];
  add_rhs = [add_rhs; kernel_coeffs.' * rhs(dofs_on_vertex)];
  add_rb = [add_rb; kernel_coeffs.' * stiff_mat(dofs_on_vertex, drchlt_dofs) * u_drchlt];

  for bv_1 = 1 : numel(add_int_dofs)
    s_vec(bv, bv_1) = kernel_coeffs.' * stiff_mat(dofs_on_vertex, space.dofs_on_vertex{add_int_dofs(bv_1).vertex_number}) * add_int_dofs(bv_1).kernel_coeffs;
  end         
end

%Stiffness matrix (including additional "interior" dofs)
A = [stiff_mat(int_dofs, int_dofs)   add_mat'; ...
                add_mat                s_vec     ];

%Right-hand side (including additional "interior" dofs)...
r = [rhs(int_dofs); ...
        add_rhs        ];
    
%...and its correction    
rb = [    stiff_mat(int_dofs, drchlt_dofs)*u_drchlt; ...
                            add_rb                     ];

%"interior" dofs = "original" interior dofs + the discared "boundary" vertex dofs                      
int_dofs_new = [int_dofs add_dofs];

%Solving the system
u_newBasis(int_dofs_new) = A \ (r - rb);
u_newBasis(drchlt_dofs) = u_drchlt;

%Switching to the usual basis
u = B_change * u_newBasis';

end
