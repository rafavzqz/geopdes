% SP_GET_LOCAL_INTERIOR_FUNCTIONS: Compute the local indices of interior basis functions on a patch.
%
% [indices] = sp_get_local_interior_functions (space, patch)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch_C1)
%    patch:   index of the patch
%
% OUTPUT:
%    indices: local indices of the C^1 interior basis functions.
%
% Copyright (C) 2022 Rafael Vazquez
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

function [interior_functions_per_patch] = sp_get_local_interior_functions (space, iptc_ind)

  [ii,jj] = ndgrid (3:space.sp_patch{iptc_ind}.ndof_dir(1)-2, 3:space.sp_patch{iptc_ind}.ndof_dir(2)-2);
  interior_functions_per_patch = sub2ind (space.sp_patch{iptc_ind}.ndof_dir, ii(:)', jj(:)');

end
