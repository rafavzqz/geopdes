% MP_SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection, in multipatch geometries.
%
%   [u, dofs] = mp_sp_drchlt_l2_proj (sp, msh, h, gnum, boundaries, refs)
%
% INPUT:
%
%  sp:         object representing the space of trial functions (see sp_bspline_2d)
%  msh:        object containing the domain partition and the quadrature rule (see msh_2d)
%  h:          function handle to compute the Dirichlet condition
%  gnum:       global numbering of the degrees of freedom (see mp_interface_2d)
%  boundaries: array of strcutures containing the information for the boundaries (see mp_geo_load)
%  refs:       boundary references on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011 Rafael Vazquez
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

function [u, dofs] = mp_sp_drchlt_l2_proj (sp, msh, h, gnum, boundaries, refs)

  dofs = [];
  ndof = max ([gnum{:}]);
  M    = spalloc (ndof, ndof, ndof);
  rhs  = spalloc (ndof, 1, ndof);

  for iref = refs
    for bnd_side = 1:boundaries(iref).nsides
      iptc = boundaries(iref).patches(bnd_side);
      iside = boundaries(iref).faces(bnd_side);

      msh_bnd = msh_eval_boundary_side (msh{iptc}, iside);
      sp_bnd  = sp_eval_boundary_side (sp{iptc}, msh_bnd);

      global_dofs = gnum{iptc}(sp_bnd.dofs);
      dofs = union (dofs, global_dofs);

      if (size (msh_bnd.geo_map, 1) == 2)
        [x, y] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                       squeeze (msh_bnd.geo_map(2,:,:)));

        hval = reshape (h (x, y, iref), sp_bnd.ncomp, msh_bnd.nqn, msh_bnd.nel);

      elseif (size (msh_bnd.geo_map, 1) == 3)
        [x, y, z] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                          squeeze (msh_bnd.geo_map(2,:,:)), ...
                          squeeze (msh_bnd.geo_map(3,:,:)));

        hval = reshape (h (x, y, z, iref), sp_bnd.ncomp, msh_bnd.nqn, msh_bnd.nel);
      end

      M_side = op_u_v (sp_bnd, sp_bnd, msh_bnd, ones (msh_bnd.nqn, msh_bnd.nel));
      M(global_dofs, global_dofs) = M(global_dofs, global_dofs) + M_side;

      rhs_side = op_f_v (sp_bnd, msh_bnd, hval);
      rhs(global_dofs) = rhs(global_dofs) + rhs_side;
    end
  end

  dofs = unique (dofs);
  u = M(dofs, dofs) \ full (rhs(dofs));

end
