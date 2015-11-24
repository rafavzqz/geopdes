% GET_BOUNDARY_INDICES: given a set of indices (elements or dofs) in the
%  whole domain, collect those that belong to a given boundary.
%
%   bnd_indices = get_boundary_indices (iside, size_dir, indices)
%
% INPUT:
%
%   iside:    number of the side where to compute the indices
%   size_dir: the size in each direction (usually given by ndof_dir or nel_dir)
%   indices:  set of indices in the whole domain
%
% OUTPUT:
%
%   bnd_indices: set of indices that belong to the boundary, with the
%     local numbering of the boundary.
%
% Copyright (C) 2015 Rafael Vazquez
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

function bnd_indices = get_boundary_indices (iside, size_dir, indices)

  ndim = numel (size_dir);
  ind2 = ceil (iside/2);
  ind = setdiff (1:ndim, ind2);

  if (mod(iside,2) == 1)
    boundary_ind = 1;
  else
    boundary_ind = size_dir (ind2);
  end

  indsub = cell (1, ndim);
  [indsub{:}] = ind2sub (size_dir, indices);
  aux = find (indsub{ind2} == boundary_ind);
  ppp = cellfun (@(x) x(aux), indsub(ind), 'UniformOutput', false);
  bnd_indices = sub2ind ([size_dir(ind), 1], ppp{:});

end
