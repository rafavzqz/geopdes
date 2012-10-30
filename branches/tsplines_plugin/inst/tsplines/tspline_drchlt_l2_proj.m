% TSPLINE_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = tspline_drchlt_l2_proj (tspline, h, sides)
%
% INPUT:
%
%  tspline: Bezier extraction read from a file of the T-splines plug-in for Rhino
%                (see read_bezier_extraction)
%  h:     function handle to compute the Dirichlet condition
%  sides: boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom, with global numbering
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011, 2012 Rafael Vazquez
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

function [u, dofs] = tspline_drchlt_l2_proj (tspline, h, sides)

  rhs  = zeros (tspline.ndof, 1);

  dofs_node_sets = [];
  dofs_side_sets = [];
  ncounter = 0;

  u = zeros (tspline.ndof, 1);
  u_node_sets = zeros (tspline.ndof, 1);
  
  for iside = sides
    if (strcmp (tspline.boundary{iside}.type, 'node'))
% For node boundaries, only constant boundary conditions are supported
      dofs_node_sets = union (dofs_node_sets, tspline.boundary{iside}.dofs);
      u_node_sets(tspline.boundary{iside}.dofs) = h (0, 0, iside);

    elseif (strcmp (tspline.boundary{iside}.type, 'side'))
% For side boundaries, we apply the L2 projection
      [msh_bnd, sp_bnd] = tspline_boundary (tspline, iside);
      dofs_side_sets = union (dofs_side_sets, sp_bnd.connectivity(:)');

      x = reshape (msh_bnd.geo_map(1,:,:), msh_bnd.nqn, msh_bnd.nel);
      y = reshape (msh_bnd.geo_map(2,:,:), msh_bnd.nqn, msh_bnd.nel);

      hval = reshape (h (x, y, iside), msh_bnd.nqn, msh_bnd.nel);

% Mass matrix to compute the L2-projection
      [rs, cs, vs] = ...
             op_u_v (sp_bnd, sp_bnd, msh_bnd, ones (msh_bnd.nqn, msh_bnd.nel));

      rows(ncounter+(1:numel(rs))) = rs;
      cols(ncounter+(1:numel(rs))) = cs;
      vals(ncounter+(1:numel(rs))) = vs;

      ncounter = ncounter + numel (rs);

% Right-hand side to compute the L2-projection
      rhs = rhs + op_f_v (sp_bnd, msh_bnd, hval);
    end
  end

% Write the two kind of boundaries (side and nodes) into one single vector
  if (ncounter ~= 0)
    M = sparse (rows, cols, vals);
    u(dofs_side_sets) = M(dofs_side_sets, dofs_side_sets) \ rhs(dofs_side_sets, 1);
  end
  u(dofs_node_sets) = u_node_sets (dofs_node_sets);

  dofs = union (dofs_node_sets, dofs_side_sets);

end
