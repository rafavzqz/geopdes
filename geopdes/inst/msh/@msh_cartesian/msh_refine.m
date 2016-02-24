% MSH_REFINE: construct a refined mesh from a given one. The function only
%                refines the mesh, it does not refine any space.
%
%     msh_fine = msh_refine (msh, nsub);
%
% INPUTS:
%     
%     msh:  an object of the msh_cartesian class (see msh_cartesian)
%     nsub: number of uniform subdivisions to apply on the elements, and for each direction
%   
% OUTPUT:
%
%     msh_fine: an object of the class msh_cartesian (see msh_cartesian)
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

function msh_fine = msh_refine (msh, nsub)

  rule = msh_gauss_nodes (msh.nqn_dir);
  [~, zeta] = kntrefine (msh.breaks, nsub-1, ones(1, msh.ndim), zeros(1, msh.ndim));
  [qn, qw] = msh_set_quad_nodes (zeta, rule);

  auxiliary_geometry.rdim = msh.rdim;
  auxiliary_geometry.map = msh.map;
  auxiliary_geometry.map_der = msh.map_der;
  if (isfield (struct (msh), 'map_der2'))
    auxiliary_geometry.map_der2 = msh.map_der2;
  end

  boundary = ~isempty (msh.boundary);
  if (msh.ndim > 1)
    bnd = [];
    for ii = 1:numel (msh.boundary)
      bnd(ii).rdim = msh.boundary(ii).rdim;
      bnd(ii).map = msh.boundary(ii).map;
      bnd(ii).map_der = msh.boundary(ii).map_der;
    end
    auxiliary_geometry.boundary = bnd;
  end
  msh_fine = msh_cartesian (zeta, qn, qw, auxiliary_geometry, 'boundary', boundary);

end
