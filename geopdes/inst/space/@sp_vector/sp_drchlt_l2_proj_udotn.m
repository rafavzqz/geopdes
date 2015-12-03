% SP_DRCHLT_L2_PROJ_UDOTN: assign the normal degrees of freedom through an L2 projection. To be used with the 'RT' spaces.
%                           The imposed condition reads   u*n = h*n
%
%   [vel, normal_dofs] = sp_drchlt_l2_proj_udotn (space, msh, bnd_sides, bnd_func)
%
% INPUTS:
%     
%    space:     space object (see sp_vector)
%    msh:       mesh object (see msh_cartesian)
%    bnd_sides: boundary sides on which the Dirichlet condition is imposed
%    bnd_func:  the condition to be imposed (h in the equations)
%   
% OUTPUT:
%
%     vel:         assigned value to the normal degrees of freedom
%     normal_dofs: global numbering of the normal basis functions
%
% Copyright (C) 2014, 2015 Rafael Vazquez
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

function [u, dofs] = sp_drchlt_l2_proj_udotn (space, msh, sides, bnd_func)

  if (~strcmpi (space.transform, 'div-preserving'))
    error ('The function only works with div-conforming spaces')
  end

  dofs = [];
  nent = 0;
  for iside = sides
    nent = nent + msh.boundary(iside).nel * space.boundary(iside).nsh_max^2;
    dofs = union (dofs, space.boundary(iside).dofs);
  end

  rows = zeros (nent, 1);
  cols = zeros (nent, 1);
  vals = zeros (nent, 1);
  rhs = zeros (space.ndof, 1);

  ncounter = 0;
  for iside = sides
    msh_side = msh_eval_boundary_side (msh, iside);
    sp_side = sp_eval_boundary_side (space, msh_side);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    g = bnd_func (x{:}, iside);

    [rs, cs, vs] = op_u_v (sp_side, sp_side, msh_side, ones(size(x{1})));

    bnd_dofs = space.boundary(iside).dofs;
    rows(ncounter+(1:numel(rs))) = bnd_dofs(rs);
    cols(ncounter+(1:numel(rs))) = bnd_dofs(cs);
    vals(ncounter+(1:numel(rs))) = vs;
    ncounter = ncounter + numel (rs);

    rhs(bnd_dofs) = rhs(bnd_dofs) + op_fdotn_v (sp_side, msh_side, g);
  end

  M = sparse (rows(1:ncounter), cols(1:ncounter), vals(1:ncounter));
  u = M(dofs, dofs) \ rhs(dofs, 1);

end
