% SP_BILAPLACIAN_DRCHLT_C1: assign the degrees of freedom of essential boundary conditions (value and normal derivative) through a projection.
%  On boundary vertices, the kernel is computed to remove linear dependencies when restricting the functions to the boundary.
%
%   [u, dofs, kernel_info] = sp_bilaplacian_drchlt_C1 (sp, msh, refs, h, dudn)
%
% INPUT:
%
%  sp:         object representing the multipatch space of trial functions (see sp_multipatch_C1)
%  msh:        object containing the domain partition and the quadrature rule (see msh_multipatch)
%  refs:       boundary references on which the conditions are imposed
%  h:          function handle to compute the Dirichlet condition
%  dudn:       function handle to compute the Neumann condition
%
% OUTPUT:
%
%  u:           assigned value to the degrees of freedom
%  dofs:        global numbering of the corresponding basis functions
%  kernel_info: a struct with information of kernel computation, containing:
%              - vertices_numbers: vertices which contain a function in the kernel
%              - all_vertex_dofs:  all functions on those vertices
%              - quasi_interior_dofs: functions that will be treated as
%                                     internal ones (as many as in the kernel)
%              - B_change_local: coefficients of the functions in the kernel,
%                                in terms of vertex basis functions. Matrix of size
%                                numel(all_vertex_dofs) x numel (quasi_interior_dofs)
%
% Copyright (C) 2022-2023 Cesare Bracco, Andrea Farahat, Rafael Vazquez
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

function [u_drchlt, drchlt_dofs, kernel_info] = sp_bilaplacian_drchlt_C1 (space, msh, refs, h, dudn)

M = spalloc (space.ndof, space.ndof, space.ndof);
rhs = zeros (space.ndof, 1);

M2 = spalloc (space.ndof, space.ndof, space.ndof);
rhs2 = zeros (space.ndof, 1);

drchlt_dofs = [];
drchlt_dofs2 = [];

boundaries = msh.boundaries;
for iref = refs
  href = @(varargin) h(varargin{:}, iref);
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);

    msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside);
    msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside);
    sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
    sp_bnd_struct = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);

    [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
    [~,icol] = find (Cpatch(sp_bnd.boundary(iside).dofs,:));
    [~,jcol] = find (Cpatch(sp_bnd.boundary(iside).adjacent_dofs,:));
    
    drchlt_dofs = union (drchlt_dofs, Cpatch_cols(icol));
    drchlt_dofs2 = union (drchlt_dofs2, Cpatch_cols(jcol));
    
    for idim = 1:msh.rdim
      x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    coeff_at_qnodes = ones (size(x{1}));
    dudn_at_qnodes = dudn (x{:},iref) .* msh_side.charlen;

    M(Cpatch_cols,Cpatch_cols) = M(Cpatch_cols, Cpatch_cols) + ...
      Cpatch.' * op_u_v (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes) * Cpatch;
    rhs(Cpatch_cols) = rhs(Cpatch_cols) + Cpatch.' * op_f_v (sp_bnd_struct, msh_side, href(x{:}));
    
    M2(Cpatch_cols,Cpatch_cols) = M2(Cpatch_cols,Cpatch_cols) + ...
      Cpatch.' * op_gradu_n_gradv_n (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes.*msh_side.charlen) * Cpatch;
    rhs2(Cpatch_cols) = rhs2(Cpatch_cols) + Cpatch.' * op_gradv_n_f (sp_bnd_struct, msh_side, dudn_at_qnodes); % I am missing the other part of the vector. It is in M2 :-)
  end
end

M_bdry = M + M2;
dofs_to_remove = [];
vertices_numbers = [];
row_indices = [];
count_vert = 0;

% Check the kernel of vertex functions on Dirichlet boundary vertices
% Pick up the basis function with the max. abs. coeff in the kernel, 
%  remove it from drchlt_dofs, and add the function in the kernel into the
%  internal part (it goes in the output)
B_change_local = [];
n_boundaries = numel(msh.boundaries); % number of boundary edges
global_refs = numel(space.interfaces) - n_boundaries + refs; % global numbering of Dirichlet boundary edges

for iv = 1 : numel(space.vertices)
  if (space.vertices(iv).boundary_vertex && space.vertices(iv).valence_p>1)
  % Loop just over Dirichlet boundary vertices
    if (~isempty(intersect(global_refs, space.vertices(iv).edges)))
      patches = space.vertices(iv).patches([1 end]);

      operations = space.vertices(iv).patch_reorientation([1 end], :);
      indices_loc_R = indices_reorientation(space.sp_patch{patches(1)}.ndof_dir, operations(1, :));
      indices_loc_L = indices_reorientation(space.sp_patch{patches(2)}.ndof_dir, operations(2, :));

      indices_loc_R = indices_loc_R(:);
      indices_loc_L = indices_loc_L(:);

%       Cpatch_ind_R = indices_loc_R([2 3 space.sp_patch{patches(1)}.ndof_dir(1)+2]);
%       Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(2)}.ndof_dir(1)+1 space.sp_patch{patches(2)}.ndof_dir(1)+2 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
      Cpatch_ind_R = indices_loc_R([1 2 3 space.sp_patch{patches(1)}.ndof_dir(1)+[1 2]]);
      if (space.vertices(iv).valence_p == 2)
        Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
      else
        Cpatch_ind_L = indices_loc_L([2 space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
      end

% TODO: it is probably enough to use CC_vertices (and CC_edges)
      [Cpatch1, Cpatch_cols1] = sp_compute_Cpatch (space, patches(1));
      [Cpatch2, Cpatch_cols2] = sp_compute_Cpatch (space, patches(2));

      [~,~,inds1] = intersect (space.dofs_on_vertex{iv}, Cpatch_cols1);
      [~,~,inds2] = intersect (space.dofs_on_vertex{iv}, Cpatch_cols2);

      M_ker = [Cpatch1(Cpatch_ind_R, inds1); ...
               Cpatch2(Cpatch_ind_L, inds2)];
%       M_ker = M_bdry(space.dofs_on_vertex{iv}, space.dofs_on_vertex{iv});
      ker = null(full(M_ker));
      if (~isempty(ker))
        count_vert = count_vert + 1;
        [~, ind] = max(abs(ker));

        row_inds = (count_vert-1)*6 + (1:6);
        B_change_local = blkdiag (B_change_local, ker(:));

        dofs_on_vertex = space.dofs_on_vertex{iv};
        vertices_numbers(count_vert) = iv;
        dofs_to_remove(count_vert) = dofs_on_vertex(ind);
        row_indices(row_inds) = dofs_on_vertex;
      end
    end
  end
end

kernel_info = struct ('vertices_numbers', vertices_numbers, 'all_vertex_dofs', row_indices, 'quasi_interior_dofs', dofs_to_remove, 'B_change_local', sparse(B_change_local));

drchlt_dofs = union (drchlt_dofs, drchlt_dofs2);
drchlt_dofs = setdiff(drchlt_dofs, dofs_to_remove);

u_drchlt = M_bdry(drchlt_dofs,drchlt_dofs) \ (rhs(drchlt_dofs) + rhs2(drchlt_dofs));

end

function indices = indices_reorientation (ndof_dir, operations)
  ndof = prod (ndof_dir);
  indices = reshape (1:ndof, ndof_dir);
  if (operations(1))
    indices = flipud (indices);
  end
  if (operations(2))
    indices = fliplr (indices);
  end
  if (operations(3))
    indices = indices.';
  end   
end