% SP_DRCHLT_L2_PROJ_UXN_3D: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%                           The imposed condition reads   uxn = h
%
%   [u, dofs] = sp_drchlt_l2_proj_uxn_3d (sp, msh, h, sides)
%
% INPUT:
%
%  sp:    structure representing the space of trial functions (see sp_bspline_2d_phys)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%  h:     function handle to compute the Dirichlet condition
%  sides: boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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


function [u, dofs] = sp_drchlt_l2_proj_uxn_3d (sp, msh, h, sides)

  dofs = unique ([sp.boundary(sides).dofs]);
  M    = spalloc (sp.ndof, sp.ndof, sp.ndof);
  rhs  = spalloc (sp.ndof, 1, numel (dofs));

  for iside = sides

    msh_bnd = msh.boundary(iside);
    sp_bnd  = sp.boundary(iside);

    [x, y, z] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                      squeeze (msh_bnd.geo_map(2,:,:)), ...
                      squeeze (msh_bnd.geo_map(3,:,:)));

    hval = reshape (h (x, y, z, iside), 3, msh_bnd.nqn, msh_bnd.nel);

    M_s = op_uxn_vxn_3d (sp_bnd, sp_bnd, msh_bnd, ones(msh_bnd.nqn, msh_bnd.nel));
    M(sp_bnd.dofs, sp_bnd.dofs) = M(sp_bnd.dofs, sp_bnd.dofs) + M_s;

    rhs_s = op_f_vxn_3d (sp_bnd, msh_bnd, hval);
    rhs(sp_bnd.dofs) = rhs(sp_bnd.dofs) + rhs_s;

  end

  u = M(dofs, dofs) \ full (rhs(dofs));

end

