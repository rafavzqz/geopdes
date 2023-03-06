% MSH_REFINE: construct a refined mesh from a given one. The function only
%                refines the mesh, it does not refine any space.
%
%     msh_fine = msh_refine (msh, nsub);
%
% INPUTS:
%     
%     msh:  an object of the msh_multipatch class (see msh_multipatch)
%     nsub: number of uniform subdivisions to apply on the elements, and for each direction
%   
% OUTPUT:
%
%     msh_fine: the refined mesh, an object of the class msh_multipatch (see msh_multipatch)
%
% Copyright (C) 2015 Rafael Vazquez
% Copyright (C) 2023 Pablo Antolin
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

  boundary = ~isempty (msh.boundary);
  msh_patch = cell (msh.npatch, 1);
  for iptc = 1:msh.npatch
    msh_ptc = msh.msh_patch{iptc};
    rule = msh_gauss_nodes (msh_ptc.nqn_dir);
    [~, zeta] = kntrefine (msh_ptc.breaks, nsub-1, ones(1, msh_ptc.ndim), zeros(1, msh_ptc.ndim));
    [qn, qw] = msh_set_quad_nodes (zeta, rule);

    auxiliary_geometry.rdim = msh_ptc.rdim;
    auxiliary_geometry.map = msh_ptc.map;
    auxiliary_geometry.map_der = msh_ptc.map_der;
    if (isfield (struct (msh_ptc), 'map_der2'))
      auxiliary_geometry.map_der2 = msh_ptc.map_der2;
    end
    if (isfield (struct (msh_ptc), 'map_der3'))
      auxiliary_geometry.map_der3 = msh_ptc.map_der3;
    end
    if (isfield (struct (msh_ptc), 'map_der4'))
      auxiliary_geometry.map_der4 = msh_ptc.map_der4;
    end
    der2 = msh.msh_patch{1}.der2;
    der3 = msh.msh_patch{1}.der3;
    der4 = msh.msh_patch{1}.der4;

    bnd = [];
    for ii = 1:numel (msh_ptc.boundary)
      bnd(ii).rdim = msh_ptc.boundary(ii).rdim;
      bnd(ii).map = msh_ptc.boundary(ii).map;
      bnd(ii).map_der = msh_ptc.boundary(ii).map_der;
      bnd(ii).map_der2 = msh_ptc.boundary(ii).map_der2;
      bnd(ii).map_der3 = msh_ptc.boundary(ii).map_der4;
      bnd(ii).map_der4 = msh_ptc.boundary(ii).map_der4;
    end
    auxiliary_geometry.boundary = bnd;
    msh_patch{iptc} = msh_cartesian (zeta, qn, qw, auxiliary_geometry, 'boundary', boundary, ...
    'der2', der2, 'der3', der3, 'der4', der4);
  end

  msh_fine = msh_multipatch (msh_patch, msh.boundaries);

end
