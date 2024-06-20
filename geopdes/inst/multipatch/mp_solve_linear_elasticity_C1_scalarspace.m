% MP_SOLVE_LINEAR_ELASTICITY_C1_SCALARSPACE: Solve a linear elasticity problem in a multipatch domain.
%  It uses an object space of type sp_multipatch_C1 with scalar values, instead of vector valued.
%
% Example to solve the linear elasticity problem
%
%      - div (sigma(u)) = f    in Omega
%      sigma(u) \cdot n = g    on Gamma_N
%                     u = h    on Gamma_D
%
% with   sigma(u) = mu*(grad(u) + grad(u)^t) + lambda*div(u)*I,
% and \Omega is an analysis-suitable G1 multipatch domain.
%
%   u:          displacement vector
%   sigma:      Cauchy stress tensor
%   lambda, mu: Lame' parameters
%   I:          identity tensor
%
% USAGE:
%
%  [geometry, msh, space, u] = 
%             mp_solve_linear_elasticity_C1_scalarspace (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (NOT IMPLEMENTED)
%    - drchlt_sides: sides with Dirichlet boundary condition (ONLY HOMOGENEOUS CONDITIONS)
%    - lambda_lame:  first Lame' parameter
%    - mu_lame:      second Lame' parameter
%    - f:            source term
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
% Copyright (C) 2010, 2011 Carlo de Falco
% Copyright (C) 2010, 2011, 2015, 2017, 2022 Rafael Vazquez
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

function [geometry, msh, space, u] = ...
              mp_solve_linear_elasticity_C1_scalarspace (problem_data, method_data)

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
mat = op_su_ev_mp (space, space, msh, lambda_lame, mu_lame);
rhs = op_f_v_mp_vector (space, msh, f);

% % % Apply Neumann boundary conditions
% % Nbnd = cumsum ([0, boundaries.nsides]);
% % for iref = nmnn_sides
% %   iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
% %   gref = @(varargin) g(varargin{:},iref);
% %   rhs_nmnn = op_f_v_mp (space.boundary, msh.boundary, gref, iref_patch_list);
% %   rhs(space.boundary.dofs) = rhs(space.boundary.dofs) + rhs_nmnn;
% % end
% % 
% % if (exist ('press_sides', 'var'))
% % for iref = press_sides
% %   rhs_press = zeros (space.boundary.ndof, 1);
% %   for iside = 1:numel(boundaries(iref).nsides)
% %     patch = boundaries(iref).patches(iside);
% %     side = boundaries(iref).faces(iside);
% %     msh_side = msh_eval_boundary_side (msh.msh_patch{patch}, side);
% %     sp_side  = sp_eval_boundary_side (space.sp_patch{patch}, msh_side);
% % 
% %     x = cell (msh_side.rdim, 1);
% %     for idim = 1:msh_side.rdim
% %       x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
% %     end
% %     pval = reshape (p (x{:}, iside), msh_side.nqn, msh_side.nel);
% % 
% %     rhs_press(space.boundary.gnum{patch}) = rhs_press(space.boundary.gnum{patch}) + op_pn_v (sp_side, msh_side, pval);
% %   end
% %   rhs(space.boundary.dofs) = rhs(space.boundary.dofs) - rhs_press;
% % end
% % end
% % 
% % symm_dofs = [];
% % if (exist ('symm_sides', 'var'))
% % for iref = symm_sides
% %   if (~strcmpi (space.transform, 'grad-preserving'))
% %     error ('The symmetry condition is only implemented for spaces with grad-preserving transform')
% %   end
% %   for iside = 1:boundaries(iref).nsides
% %     patch = boundaries(iref).patches(iside);
% %     side = boundaries(iref).faces(iside);
% %     msh_side = msh_eval_boundary_side (msh.msh_patch{patch}, side);
% %     normal_comp = zeros (msh_side.rdim, msh_side.nqn * msh_side.nel);
% %     for idim = 1:msh_side.rdim
% %       normal_comp(idim,:) = reshape (msh_side.normal(idim,:,:), 1, msh_side.nqn*msh_side.nel);
% %     end
% %     
% %     parallel_to_axes = false;
% %     for ind = 1:msh_side.rdim
% %       ind2 = setdiff (1:msh_side.rdim, ind);
% %       if (all (all (abs (normal_comp(ind2,:)) < 1e-10)))
% %         bnd_side = space.sp_patch{patch}.boundary(side);
% %         dofs = bnd_side.dofs(bnd_side.comp_dofs{ind});
% %         symm_dofs = union (symm_dofs, space.gnum{patch}(dofs));
% %         parallel_to_axes = true;
% %         break
% %       end
% %     end
% %     if (~parallel_to_axes)
% %       error ('mp_solve_linear_elasticity: We have only implemented the symmetry condition for boundaries parallel to the axes')
% %     end
% %   end
% % end
% % end

% Apply Dirichlet boundary conditions
% TODO: implement non-homogeneous conditions
drchlt_dofs = [];
for iref = 1:numel(drchlt_sides)
  bnd_ref = drchlt_sides(iref);
  scalar_dofs_on_ref = [];
  for bnd_side = 1:msh.boundaries(bnd_ref).nsides
    iptc = msh.boundaries(bnd_ref).patches(bnd_side);
    iside = msh.boundaries(bnd_ref).faces(bnd_side);

%     msh_side = msh.msh_patch{iptc}.boundary(iside);
    side_dofs = space.sp_patch{iptc}.boundary(iside).dofs;

    [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);

    [~,scalar_dofs] = find (Cpatch(side_dofs,:));
    scalar_dofs_on_ref = union (scalar_dofs_on_ref, Cpatch_cols(scalar_dofs));
  end
  for icomp = 1:msh.rdim
    drchlt_dofs = union (drchlt_dofs, (icomp-1)*space.ndof + scalar_dofs_on_ref);
  end
end

ndof = msh.rdim * space.ndof;
u = zeros (ndof, 1);

int_dofs = setdiff (1:ndof, drchlt_dofs);
% Solve the linear system
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

end
