% MSH_MULTIPATCH: constructor of the msh class for multipatch grids.
%
%     msh = msh_multipatch (meshes, geometry);
%     msh = msh_multipatch (meshes, boundaries, interfaces);
%
% INPUTS:
%     
%     meshes:     cell-array of Cartesian meshes (see msh_cartesian)
%     interfaces: structure with the information of the interfaces between patches (see mp_geo_load)
%     boundaries:
%   
% OUTPUT:
%
%     msh: object containing the following fields and methods
%
%     FIELD_NAME    (SIZE)                    DESCRIPTION
%     npatch        (scalar)                  the number of patches in the domain
%     ndim          (scalar)                  dimension of the parametric space
%     rdim          (scalar)                  dimension of the physical space
%     nel           (scalar)                  total number of elements of the partition
%     nel_per_patch (1 x npatch array)        number of elements on each patch
%     nqn           (scalar)                  number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)         number of quadrature nodes per element in each parametric direction
%     msh_patch     (1 x npatch)
%     boundary      () it contains an (ndim-1)-dimensional 'msh_cartesian' object for each side of the boundary (only when boundary is set to true)
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

function msh = msh_multipatch (meshes)

  aux = struct ([meshes{:}]);
  
  msh.npatch = numel (meshes);
  msh.ndim = meshes{1}.ndim;
  msh.rdim = meshes{1}.rdim;
  
  msh.nel = sum ([aux.nel]);
  msh.nel_per_patch = [aux.nel];
  msh.msh_patch = meshes;
  
% All the boundary info is still missing  
  msh = class (msh, 'msh_multipatch');  
end
