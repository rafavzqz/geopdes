% SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = sp_drchlt_l2_proj (sp, msh, h, sides)
%
% INPUT:
%
%  sp:    object defining the space of discrete functions (see sp_bspline_2d)
%  msh:   object defining the domain partition and the quadrature rule (see msh_2d)
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

function [u, dofs] = sp_drchlt_l2_proj (sp, msh, h, sides)

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
  for iside = sides
    
    msh_bnd = msh_eval_boundary_side (msh, iside);
    sp_bnd  = sp_eval_boundary_side (sp, msh_bnd);

    if (size (msh_bnd.geo_map, 1) == 2)
      [x, y] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                     squeeze (msh_bnd.geo_map(2,:,:)));

      hval = reshape (h (x, y, iside), sp.ncomp, msh_bnd.nqn, msh_bnd.nel);

    elseif (size (msh_bnd.geo_map, 1) == 3)
      [x, y, z] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                        squeeze (msh_bnd.geo_map(2,:,:)), ...
                        squeeze (msh_bnd.geo_map(3,:,:)));

      hval = reshape (h (x, y, z, iside), sp_bnd.ncomp, msh_bnd.nqn, msh_bnd.nel);
    end

    [rs, cs, vs] = ...
             op_u_v (sp_bnd, sp_bnd, msh_bnd, ones (msh_bnd.nqn, msh_bnd.nel));

    rows(ncounter+(1:numel(rs))) = sp_bnd.dofs(rs);
    cols(ncounter+(1:numel(rs))) = sp_bnd.dofs(cs);
    vals(ncounter+(1:numel(rs))) = vs;
    ncounter = ncounter + numel (rs);

    rhs_side = op_f_v (sp_bnd, msh_bnd, hval);
    rhs(sp_bnd.dofs) = rhs(sp_bnd.dofs) + rhs_side;
  end

  M = sparse (rows, cols, vals);
  u = M(dofs, dofs) \ rhs(dofs, 1);

end
