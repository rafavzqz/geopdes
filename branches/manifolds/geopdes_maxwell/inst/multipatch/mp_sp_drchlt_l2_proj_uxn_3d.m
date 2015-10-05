% MP_SP_DRCHLT_L2_PROJ_UXN_3D: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%                           The imposed condition reads   uxn = h
%
%   [u, dofs] = mp_sp_drchlt_l2_proj_uxn_3d (sp, msh, h, gnum, ornt, boundaries, refs)
%
% INPUT:
%
%  sp:    structure representing the space of trial functions (see sp_bspline_2d_phys)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%  h:     function handle to compute the Dirichlet condition
%  gnum:  global numbering of the degrees of freedom (see mp_interface_2d)
%  ornt:  orientation of the degrees of freedom (see mp_interface_hcurl_3d)
%  boundaries: array of strcutures containing the information for the boundaries (see mp_geo_load)
%  refs:  boundary references on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2012 Rafael Vazquez
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

function [u, dofs] = mp_sp_drchlt_l2_proj_uxn_3d (sp, msh, h, gnum, ornt, boundaries, refs)

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


      [x, y, z] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                        squeeze (msh_bnd.geo_map(2,:,:)), ...
                        squeeze (msh_bnd.geo_map(3,:,:)));

      hval = reshape (h (x, y, z, iref), 3, msh_bnd.nqn, msh_bnd.nel);

      M_side = op_uxn_vxn_3d (sp_bnd, sp_bnd, msh_bnd, ones(size(x)));
      ornt_bnd = ornt{iptc}(sp_bnd.dofs);
      ornt_matrix = spdiags (ornt_bnd', 0, sp_bnd.ndof, sp_bnd.ndof);
      M_side = ornt_matrix * M_side * ornt_matrix;
      M(global_dofs, global_dofs) = M(global_dofs, global_dofs) + M_side;

      rhs_side = op_f_vxn_3d (sp_bnd, msh_bnd, hval) .* ornt_bnd';
      rhs(global_dofs) = rhs(global_dofs) + rhs_side;

    end
  end

  dofs = unique (dofs);
  u = M(dofs, dofs) \ full (rhs(dofs));

end

