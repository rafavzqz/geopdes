% SOLVE_LAPLACE_COLLOCATION: Solve a Laplace problem with a B-spline discretization (non-isoparametric approach). 
%
% The function solves the diffusion problem
%
%    - div (grad (u)) = f    in Omega = F((0,1)^n)
%                   u = 0    on Gamma
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_laplace_collocation (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (not implemented yet)
%    - h:            function for Dirichlet boundary condition (not implemented yet)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - pts_case:   the choice of collocation points, 1:uniform, 2:Greville
%    - ncoll_pts:  number of collocation points in each direction, for the choice pts_case=1
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      auxiliary "mesh" object with the collocation points (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  u:        the computed degrees of freedom
%
% Since GeoPDEs is based on a mesh structure, we generate an auxiliary
%  "mesh" with only one collocation point per element, using the 
%  midpoints between collocation points as the mesh lines.
% Although this is not the most efficient way to implement collocation, it
%  allows us to use all the functionality in GeoPDEs, in particular the
%  evaluation of basis functions, without further changes.
%
% Copyright (C) 2016, 2017 Lorenzo Tamellini, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh_coll, space, u] = solve_laplace_collocation (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure, perform degree elevation and knot insertion
geometry = geo_load (geo_name);

if (isfield (geometry, 'nurbs'))
  ndim = numel (geometry.nurbs.order);
end

degelev  = max (degree - (geometry.nurbs.order-1), 0);
if (any (degelev < 0))
  warning('The degree provided is lower than the degree of the original geometry')
elseif (any (method_data.degree < 2) || any (method_data.regularity < 1))
  error ('Collocation for the Laplacian requires at least C^1 continuity. Degree should be at least 2')
end
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[knots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Compute the collocation points
coll_pts = cell (1,ndim);
switch pts_case
  case 1
    if (numel (ncoll_pts) == 1)
      ncoll_pts = ncoll_pts * ones (ndim, 1);
    end
    for idim = 1:ndim
      coll_pts{idim} = linspace (0,1,ncoll_pts);
    end
  case 2
    for idim = 1:ndim
      coll_pts{idim} = aveknt (knots{idim}, nurbs.order(idim)); 
    end
  otherwise
    error ('That choice of the collocation points is not implemented (yet)')
end

% Generate the auxiliary mesh object, with one collocation point per
%  element, and the associated space object
brk = cell (1,ndim);
for idim = 1:ndim
  coll_pts{idim} = coll_pts{idim}(:)';
  if (numel(coll_pts{idim}) > 1)
    brk{idim} = [knots{idim}(1), coll_pts{idim}(1:end-1) + diff(coll_pts{idim})/2, knots{idim}(end)];
  else
    brk{idim} = [knots{idim}(1) knots{idim}(end)];
  end
end

msh_coll = msh_cartesian (brk, coll_pts, [], geometry, 'boundary', true,'der2',true);
space = sp_nurbs (geometry.nurbs, msh_coll);

% Evaluate each spline in the collocation points and collect everything into the design matrix, i.e. a
% matrix whose n-th row is the equation collocated in the n-th point. Remember that we have built msh_coll
% to have 1 pt per element, so the number of coll pts is exactly the number of elements of msh_coll
tot_nb_coll_pts = msh_coll.nel;
nb_dofs = space.ndof;
A = sparse (tot_nb_coll_pts, nb_dofs);

% Evaluate parameterization and basis functions
coll_mesh_eval = msh_precompute (msh_coll);
sp_evals  = sp_precompute (space, coll_mesh_eval, 'gradient', false, 'laplacian', true);

% Assemble the matrix of collocation
for iel = 1:msh_coll.nel
  list_fun = sp_evals.connectivity(:,iel);
  A(iel,list_fun) = -sp_evals.shape_function_laplacians(1,:,iel);
end

% Generate the right-hand side, using the coordinates of the collocation points.
x = cell (1,ndim);
for idim = 1:ndim
  x{idim} = reshape (coll_mesh_eval.geo_map(idim,:,:), msh_coll.nel, 1);
end
rhs = problem_data.f(x{:});

% Apply homogeneous Dirichlet boundary condition
nb_boundaries = numel (sp_evals.boundary);
boundary_dofs = [];
for bb = 1:nb_boundaries
  boundary_dofs = union (boundary_dofs, sp_evals.boundary(bb).dofs);
end
internal_dofs = setdiff (1:sp_evals.ndof, boundary_dofs);
A(:,boundary_dofs) = [];

% If the matrix is not square, solve with least squares for now A'*A
if (size (A,1) ~= size (A,2))
  rhs = A' * rhs;
  A   = A' * A;
end

u = zeros (space.ndof,1);
u(internal_dofs) = A\rhs;

end
