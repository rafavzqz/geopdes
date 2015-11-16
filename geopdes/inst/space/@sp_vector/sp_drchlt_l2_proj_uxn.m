% SP_DRCHLT_L2_PROJ_UXN: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%                           The imposed condition reads   uxn = h
%
%   [u, dofs] = sp_drchlt_l2_proj_uxn (sp, msh, h, sides)
%
% INPUT:
%
%  sp:    object defining the space of discrete functions (see sp_vector)
%  msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%  h:     function handle to compute the Dirichlet condition
%  sides: boundary sides on which a Dirichlet condition is imposed
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

function [u, dofs] = sp_drchlt_l2_proj_uxn (sp, msh, h, sides)

  if (~strcmpi (sp.transform, 'curl-preserving'))
    warning ('SP_DRCHLT_L2_PROJ_UXN is deprecated. Consider using SP_DRCHLT_L2_PROJ, instead')
  end

  rhs  = zeros (sp.ndof, 1);

  dofs = [];
  nent = 0;
  for iside = sides
    nent = nent + msh.boundary(iside).nel * sp.boundary(iside).nsh_max^2;
    dofs = union (dofs, sp.boundary(iside).dofs);
  end

  rows = zeros (nent, 1);
  cols = zeros (nent, 1);
  vals = zeros (nent, 1);
  
  ncounter = 0;

  if (sp.ncomp == 2)
    for iside = sides
      msh_bnd = msh_eval_boundary_side (msh, iside);
      sp_bnd  = sp_eval_boundary_side (sp, msh_bnd);

      x = squeeze (msh_bnd.geo_map(1,:,:));
      y = squeeze (msh_bnd.geo_map(2,:,:));

      hval = reshape (h (x, y, iside), msh_bnd.nqn, msh_bnd.nel);

%       [rs, cs, vs] = op_uxn_vxn_2d (sp_bnd, sp_bnd, msh_bnd, ones(msh_bnd.nqn, msh_bnd.nel));
      [rs, cs, vs] = op_u_v (sp_bnd, sp_bnd, msh_bnd, ones(msh_bnd.nqn, msh_bnd.nel));

      rows(ncounter+(1:numel(rs))) = sp_bnd.dofs(rs);
      cols(ncounter+(1:numel(rs))) = sp_bnd.dofs(cs);
      vals(ncounter+(1:numel(rs))) = vs;
      ncounter = ncounter + numel (rs);

      rhs_side = op_f_vxn_2d (sp_bnd, msh_bnd, hval);
      rhs(sp_bnd.dofs) = rhs(sp_bnd.dofs) + rhs_side;
    end
  elseif (sp.ncomp == 3)
    for iside = sides
      msh_bnd = msh_eval_boundary_side (msh, iside);
      sp_bnd  = sp_eval_boundary_side (sp, msh_bnd);

      x = squeeze (msh_bnd.geo_map(1,:,:));
      y = squeeze (msh_bnd.geo_map(2,:,:));
      z = squeeze (msh_bnd.geo_map(3,:,:));

      hval = reshape (h (x, y, z, iside), 3, msh_bnd.nqn, msh_bnd.nel);

      [rs, cs, vs] = op_uxn_vxn_3d (sp_bnd, sp_bnd, msh_bnd, ones(msh_bnd.nqn, msh_bnd.nel));

      rows(ncounter+(1:numel(rs))) = sp_bnd.dofs(rs);
      cols(ncounter+(1:numel(rs))) = sp_bnd.dofs(cs);
      vals(ncounter+(1:numel(rs))) = vs;
      ncounter = ncounter + numel (rs);

      rhs_s = op_f_vxn_3d (sp_bnd, msh_bnd, hval);
      rhs(sp_bnd.dofs) = rhs(sp_bnd.dofs) + rhs_s;
    end
  end

  M = sparse (rows, cols, vals);
  u = M(dofs, dofs) \ full (rhs(dofs));

end
