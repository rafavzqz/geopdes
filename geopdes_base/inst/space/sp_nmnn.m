% SP_NMNN: compute the right-hand side for Neumann boundary conditions.
%
%   rhs = sp_nmnn (sp, msh, g, sides)
%
% INPUT:
%
%  sp:    structure representing the space of trial functions (see sp_bspline_2d_phys)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%  g:     function handle to compute the Neumann condition
%  sides: boundary sides on which a Neumann condition is imposed
%
% OUTPUT:
%
%  rhs:  right-hand side for the Neumann conditions
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

function rhs = sp_nmnn (sp, msh, g, nmnn_sides)

  rhs = zeros (sp.ndof, 1);

  if (size (msh.boundary(1).geo_map, 1) == 2)
    for iside = nmnn_sides
      [x, y] = deal (squeeze (msh.boundary(iside).geo_map(1,:,:)), ...
                     squeeze (msh.boundary(iside).geo_map(2,:,:)));
      gval = reshape (g (x, y, iside), msh.boundary(iside).nqn, msh.boundary(iside).nel);

      rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + ...
          op_f_v (sp.boundary(iside), msh.boundary(iside), gval);
    end


  elseif (size (msh.boundary(1).geo_map, 1) == 3)
    for iside = nmnn_sides
      [x, y, z] = deal (squeeze (msh.boundary(iside).geo_map(1,:,:)), ...
                        squeeze (msh.boundary(iside).geo_map(2,:,:)), ...
                        squeeze (msh.boundary(iside).geo_map(3,:,:)));

      gval = reshape (g (x, y, z, iside), msh.boundary(iside).nqn, msh.boundary(iside).nel);

      rhs(sp.boundary(iside).dofs) = rhs(sp.boundary(iside).dofs) + ...
          op_f_v (sp.boundary(iside), msh.boundary(iside), gval);
    end

  end

end
