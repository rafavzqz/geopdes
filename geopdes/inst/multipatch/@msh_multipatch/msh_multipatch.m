% MSH_MULTIPATCH: constructor of the msh class for multipatch grids.
%
%     msh = msh_multipatch (meshes);
%     msh = msh_multipatch (meshes, boundaries);
%
% INPUTS:
%     
%     meshes:     cell-array of Cartesian meshes (see msh_cartesian)
%     boundaries: array of structures with the patch and side number for
%                   the domain boundaries (see mp_geo_load)
%   
% OUTPUT:
%
%     msh: object containing the following fields and methods
%
%     FIELD_NAME    (SIZE)                    DESCRIPTION
%     npatch        (scalar)                   the number of patches in the domain
%     ndim          (scalar)                   dimension of the parametric space
%     rdim          (scalar)                   dimension of the physical space
%     nel           (scalar)                   total number of elements of the partition
%     nel_per_patch (1 x npatch array)         number of elements on each patch
%     nqn           (scalar)                   number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)          number of quadrature nodes per element in each parametric direction
%     msh_patch     (1 x npatch cell-array)    the input meshes, a mesh object for each patch (see msh_cartesian)
%     boundary      (1 x 1 object)             it contains an (ndim-1)-dimensional 'msh_multipatch' object for the whole boundary
%     patch_numbers (1 x npatch array)         only for boundary objects, the volumetric patch to which the boundary patch belongs
%     side_numbers  (1 x npatch array)         only for boundary objects, the side number that the patch occupies in the volumetric patch
%     boundaries    (struct array)             information (patch and side) for each group of boundary sides (see mp_geo_load)
%
%     METHOD NAME
%     msh_evaluate_element_list: compute the parameterization (and its derivatives) at
%                       the quadrature points of a given list of elements.
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

function msh = msh_multipatch (meshes, boundaries)

  aux = struct ([meshes{:}]);
  
  msh.npatch = numel (meshes);
  msh.ndim = meshes{1}.ndim;
  msh.rdim = meshes{1}.rdim;
  
  msh.nel = sum ([aux.nel]);
  msh.nel_per_patch = [aux.nel];
  msh.msh_patch = meshes;

  msh.boundaries = [];
  msh.patch_numbers = [];
  msh.side_numbers  = [];
  msh.boundary = [];
  
  if (nargin == 2)
    msh.boundaries = boundaries;
    if (~isempty (meshes{1}.boundary))

      patch_numbers = (vertcat (boundaries.patches)).';
      side_numbers  = (vertcat (boundaries.faces)).';
      for ind = 1:numel(patch_numbers)
        msh_bnd{ind} = meshes{patch_numbers(ind)}.boundary(side_numbers(ind));
      end
      msh.boundary = msh_multipatch (msh_bnd);
      msh.boundary.patch_numbers = patch_numbers;
      msh.boundary.side_numbers  = side_numbers;
    end
  end

  msh = class (msh, 'msh_multipatch');
  
end
