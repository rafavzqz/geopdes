% SP_COMPUTE_CPATCH: Compute the matrix for B-spline representation, and
%  the indices of C^1 functions that do not vanish on a given patch.
%
% [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, patch)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch_C1)
%    patch:   index of the patch
%
% OUTPUT:
%    Cpatch:      coefficients of the linear combination of basis functions as standard B-splines,
%                   The matrix has size ndof_per_patch(iptc) x numel(Cpatch_cols).
%    Cpatch_cols: indices of the C^1 basis functions that do not vanish on the patch
%                   (see also sp_get_functions_on_patch)
%
% Copyright (C) 2022 Cesare Bracco, Rafael Vazquez
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

function [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc_ind)

  numel_interior_dofs = space.ndof_interior_per_patch;
  ndof_per_interface = cellfun (@numel, space.dofs_on_edge);
  ndof_per_vertex = cellfun (@numel, space.dofs_on_vertex);

% Interior basis functions
  global_indices = sum (numel_interior_dofs(1:iptc_ind-1)) + (1:numel_interior_dofs(iptc_ind));
  rows = sp_get_local_interior_functions (space, iptc_ind);%space.interior_dofs_per_patch{iptc_ind}; 
  cols = 1:numel_interior_dofs(iptc_ind);%global_indices; 
  vals = ones (numel_interior_dofs(iptc_ind), 1);
  Cpatch = sparse (rows, cols, vals, space.ndof_per_patch(iptc_ind), numel_interior_dofs(iptc_ind));
  Cpatch_cols = global_indices;

% Edge basis functions
  for intrfc = 1:numel(space.interfaces)
    global_indices = space.ndof_interior + sum(ndof_per_interface(1:intrfc-1)) + (1:ndof_per_interface(intrfc));
    patches = [space.interfaces(intrfc).patch1 space.interfaces(intrfc).patch2];
    [patch_is_on_interface,local_patch] = ismember (iptc_ind, patches);
    if (patch_is_on_interface)
      indices = size(Cpatch,2) + (1:numel(global_indices));
      Cpatch(:,indices) = space.CC_edges{local_patch,intrfc};
      Cpatch_cols = union (Cpatch_cols, global_indices);
    end
  end

% Vertex basis functions
  for ivrt = 1:numel(space.vertices)
    global_indices = space.ndof_interior + space.ndof_edges + sum(ndof_per_vertex(1:ivrt-1)) + (1:ndof_per_vertex(ivrt));
    [patch_is_on_vertex,local_patch] = ismember (iptc_ind, space.vertices(ivrt).patches);
    if (patch_is_on_vertex)
      indices = size(Cpatch,2) + (1:numel(global_indices));
      Cpatch(:,indices) = space.CC_vertices{ivrt}{local_patch};
      Cpatch_cols = union (Cpatch_cols, global_indices);
    end
  end

end