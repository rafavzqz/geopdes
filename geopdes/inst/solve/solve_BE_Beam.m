% SOLVE_BE_BEAM: Solve a Bernoulli-Euler beam static problem on a one-dimensional domain.
%
%      - (E(x)I(x)w'')'' = f(x)     in Omega = F(0,1)
%
% w=0  (1st Dirichlet BC) or -E*I*w'''=P    (1st Neumann BC) on x=0, x=L
% w'=0 (2st Dirichlet BC) or -E*I*w''=M_0,L (2nd Neumann BC) on x=0, x=L  
%
%   w:    deflection
%   E:    Young's modulus E
%   I:    moment of inertia
%
% USAGE:
%
%  [geometry, msh, space, w] = solve_BE_Beam (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%      fields containing information about the boundary conditions on 
%      both beam ends in the following form: [boolean boolean]
%    - drchlt1_sides: ends with 1st type Dirichlet boundary condition (Homogeneous)
%    - drchlt2_sides: ends with 2nd type Dirichlet boundary condition (Homogeneous)
%      (at least two Dirichlet BCs must be defined and at least one of the 1st type)
%    - nmnn1_sides:   ends with 1st type Neumann boundary condition
%    - nmnn2_sides:   ends with 2nd type Neumann boundary condition
%    - EI:            E(x)*I(x) function
%    - P:             value of the concentrated force (1st Neumann boundary condition)
%    - M_L and M_0:   values of the moments (2nd Neumann boundary condition)
%    - f:             function of the distributed loading
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
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_2d)
%  space:    space object that defines the discrete basis functions (see sp_vector_2d)
%  w:        the computed degrees of freedom
%
% See also ex_BE_beam_static_1 and ex_BE_beam_static_2 for example.
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2016 Viacheslav Balobanov
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, sp, u] = solve_BE_Beam (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% BC's checking
if sum(drchlt1_ends) < 1 ||... 
   sum(drchlt1_ends) + sum(drchlt2_ends) < 2 ||...
   drchlt1_ends(1) + nmnn1_ends(1) ~= 1 ||...   
   drchlt1_ends(2) + nmnn1_ends(2) ~= 1 ||...
   drchlt2_ends(1) + nmnn2_ends(1) ~= 1 ||...
   drchlt2_ends(2) + nmnn2_ends(2) ~= 1
      error('Boundary conditions are not correct')
end

% Construct geometry structure
geometry = geo_load (geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry, 'boundary', true, 'der2', true);

% Construct space structure
sp = sp_bspline (knots, degree, msh);

% Assemble the stiffness matrix
stiff_mat = op_gradgradu_gradgradv_tp (sp, sp, msh, EI);

%Assemble the force vector
rhs = op_f_v_tp (sp, msh, f);

% Apply homogeneous 1st Dirichlet boundary conditions
%  and 1st Neumann boundary conditions
drchlt_dofs = [];
for iside = 1:2*msh.ndim
  if (drchlt1_ends(iside))
    drchlt_dofs = [drchlt_dofs, sp.boundary(iside).dofs];
  elseif (nmnn1_ends(iside))
    rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + P;
  end
end

% Apply 2nd Dirichlet boundary conditions by using the Lagrange multipliers method
%  and 2nd Neumann boundary conditions
n_d2 = sum (drchlt2_ends);
C = sparse (numel(drchlt2_ends), sp.ndof);
for iside = 1:2*msh.ndim
  if (drchlt2_ends(iside) || nmnn2_ends(iside))
    msh_aux = msh_boundary_side_from_interior (msh, iside);
    sp_side = sp_precompute (sp.constructor (msh_aux), msh_aux, 'gradient', true);
    
    if (drchlt2_ends(iside))
      C(iside, sp_side.connectivity) = reshape (sp_side.shape_function_gradients(:,:,:), 1, sp_side.nsh);
    elseif (nmnn2_ends(iside))
      rhs(sp_side.connectivity) = rhs(sp_side.connectivity) - ...
          M(iside) * reshape (sp_side.shape_function_gradients(:,:,:), sp_side.nsh, 1);
    end
  end
end
stiff_mat = [stiff_mat,         C(drchlt2_ends, :).'; ...
             C(drchlt2_ends,:), sparse(n_d2, n_d2)];
rhs(sp.ndof+(1:n_d2)) = 0;
u = zeros (sp.ndof + n_d2, 1);
int_dofs = setdiff (1:(sp.ndof+n_d2), drchlt_dofs);
    
% Solve the static problem
K = stiff_mat(int_dofs, int_dofs);
F = rhs(int_dofs);
u(int_dofs) = K\F;
u = u(1:sp.ndof);

end