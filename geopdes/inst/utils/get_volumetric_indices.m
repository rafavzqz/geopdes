% GET_VOLUMETRIC_INDICES: given a set of indices (elements or dofs) on the
%  boundary, compute the corresponding indices in the whole domain.
%
%   vol_indices = get_volumetric_indices (iside, size_dir, indices)
%
% INPUT:
%
%   iside:    number of the side where to compute the indices
%   size_dir: the size in each direction (usually given by ndof_dir or nel_dir)
%   indices:  set of indices on the boundary
%
% OUTPUT:
%
%   vol_indices: corresponding indices in the global domain
%
%   For single patch geometries, this does the same work as the "dofs" field
%    in the space object.
%
% Copyright (C) 2018 Rafael Vazquez
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

function vol_indices = get_volumetric_indices (iside, size_dir, indices)

  ndim = numel (size_dir);
  ind2 = ceil (iside/2);
  ind = setdiff (1:ndim, ind2);

  if (mod(iside,2) == 1)
    boundary_ind = 1;
  else
    boundary_ind = size_dir (ind2);
  end

  indsub = cell (1, ndim);
  [indsub{ind}] = ind2sub ([size_dir(ind), 1], indices);
  indsub{ind2} = boundary_ind * ones (size (indsub{ind(1)}));

  vol_indices = sub2ind ([size_dir, 1], indsub{:});
  
end
