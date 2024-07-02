% SP_DRCHLT_C1_SHELLS: assign the degrees of freedom to impose zero displacement boundary conditions.
%  On boundary vertices, the kernel is computed to remove linear dependencies when restricting the functions to the boundary.
%
%   [u, dofs, kernel_info] = sp_drchlt_C1_shells (sp, msh, refs, drchlt_components)
%
% INPUT:
%
%  sp:         object representing the multipatch space of trial functions (see sp_multipatch_C1)
%  msh:        object containing the domain partition and the quadrature rule (see msh_multipatch)
%  refs:       boundary references on which the conditions are imposed
%  drchlt_components: cell-array, the components that are set to zero for each drchlt_side
%
% OUTPUT:
%
%  u:           assigned value to the degrees of freedom (zero for now)
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

function [u_drchlt, drchlt_dofs, kernel_info] = sp_drchlt_C1_shells (space, msh, refs, drchlt_components)

drchlt_dofs = [];

boundaries = msh.boundaries;
for iref = 1:numel(refs)
%   href = @(varargin) h(varargin{:}, iref);
  if (~exist('drchlt_components','var'))
    components = 1:3;
  else
    components = drchlt_components{iref};
  end
  scalar_dofs_on_ref = [];
  for bnd_side = 1:boundaries(refs(iref)).nsides
    iptc = boundaries(refs(iref)).patches(bnd_side);
    iside = boundaries(refs(iref)).faces(bnd_side);

    msh_side = msh.msh_patch{iptc}.boundary(iside);
    sp_bnd = space.sp_patch{iptc}.boundary(iside);

    [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);

    [~,scalar_dofs] = find (Cpatch(sp_bnd.dofs,:));
    scalar_dofs_on_ref = union (scalar_dofs_on_ref, Cpatch_cols(scalar_dofs));
%     Cpatch_bnd = space.Cpatch{iptc}(sp_bnd.dofs,:);
%     [~,scalar_dofs] = find (abs(Cpatch_bnd)>1e-15);
%     
%     scalar_dofs_on_ref = union (scalar_dofs_on_ref, scalar_dofs);
  end
  for icomp = components
    drchlt_dofs = union (drchlt_dofs, (icomp-1)*space.ndof + scalar_dofs_on_ref);
  end
end

dofs_to_remove = [];
vertices_numbers = [];
row_indices = [];
count_vert = 0;
count_fun = 0;

% Check the kernel of vertex functions on Dirichlet boundary vertices
% Pick up the basis function with the max. abs. coeff in the kernel, 
%  remove it from drchlt_dofs, and add the function in the kernel into the
%  internal part (it goes in the output)
B_change_local = [];
n_boundaries = numel(msh.boundaries); % number of boundary edges
global_refs = numel(space.interfaces) - n_boundaries + refs; % global numbering of Dirichlet boundary edges

for iv = 1 : numel(space.vertices)
  % Loop just over Dirichlet boundary vertices
  if ~isempty(intersect(global_refs, space.vertices(iv).edges))
    if (space.vertices(iv).boundary_vertex)
      patches = space.vertices(iv).patches([1 end]);

      operations = space.vertices(iv).patch_reorientation([1 end], :);
      indices_loc_R = indices_reorientation(space.sp_patch{patches(1)}.ndof_dir, operations(1, :));
      indices_loc_L = indices_reorientation(space.sp_patch{patches(2)}.ndof_dir, operations(2, :));
      indices_loc_R = indices_loc_R(:);
      indices_loc_L = indices_loc_L(:);

      Cpatch_ind_R = indices_loc_R([1 2 3]);
      Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(1)}.ndof_dir(1)+1 2*space.sp_patch{patches(1)}.ndof_dir(1)+1]);
%       Cpatch_ind_R = indices_loc_R([1 2 3 space.sp_patch{patches(1)}.ndof_dir(1)+[1 2]]);
%       if (space.vertices(iv).valence_p == 2)
%         Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
%       else
%         Cpatch_ind_L = indices_loc_L([2 space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
%       end

      [Cpatch1, Cpatch_cols1] = sp_compute_Cpatch (space, patches(1));
      [Cpatch2, Cpatch_cols2] = sp_compute_Cpatch (space, patches(2));

      [~,~,inds1] = intersect (space.dofs_on_vertex{iv}, Cpatch_cols1);
      [~,~,inds2] = intersect (space.dofs_on_vertex{iv}, Cpatch_cols2);

      M_ker = [Cpatch1(Cpatch_ind_R, inds1); ...
               Cpatch2(Cpatch_ind_L, inds2)];
%       M_ker = [space.Cpatch{patches(1)}(Cpatch_ind_R, space.dofs_on_vertex{iv}); ...
%                space.Cpatch{patches(2)}(Cpatch_ind_L, space.dofs_on_vertex{iv})];

      ker = null(full(M_ker));
      if (~isempty(ker))
        nfun = size(ker,2);
        [~, ind] = max(abs(ker)); % TODO: NOT A GOOD CHOICE (it may be repeated)

        row_inds = count_vert*6 + (1:6);
        B_change_local = blkdiag (B_change_local, ker);

        dofs_on_vertex = space.dofs_on_vertex{iv};
        vertices_numbers(count_fun+(1:nfun)) = iv;
        dofs_to_remove(count_fun+(1:nfun)) = dofs_on_vertex(ind);
        row_indices(row_inds) = dofs_on_vertex;
        count_vert = count_vert + 1;
        count_fun  = count_fun + nfun;
      end
    end
  end
end

kernel_info = struct ('vertices_numbers', vertices_numbers, 'all_vertex_dofs', row_indices, 'quasi_interior_dofs', dofs_to_remove, 'B_change_local', sparse(B_change_local));

dofs_to_remove = [dofs_to_remove(:), dofs_to_remove(:)+space.ndof, dofs_to_remove(:)+2*space.ndof];
drchlt_dofs = setdiff(drchlt_dofs, dofs_to_remove);
% drchlt_dofs = [drchlt_dofs(:); drchlt_dofs(:)+space.ndof; drchlt_dofs(:)+2*space.ndof];

u_drchlt = zeros (numel(drchlt_dofs), 1);

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
