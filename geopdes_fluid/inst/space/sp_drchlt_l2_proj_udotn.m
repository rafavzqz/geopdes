% SP_DRCHLT_L2_PROJ_UDOTN: assign the normal degrees of freedom thtough an L2 projection. To be used with the 'RT' spaces.
%                           The imposed condition reads   u*n = h*n
%
%   [vel, normal_dofs] = sp_drchlt_l2_proj_udotn (space, msh, geometry, der2, bnd_sides, bnd_func, coeff, Cpen)
%
% INPUTS:
%     
%    space:     space object (see sp_vector_piola_transform_2d)
%    msh:       mesh object (see msh_2d)
%    geometry:  geometry object (see geo_load)
%    der2:      a logical telling whether to compute the second derivative.
%               For Raviart-Thomas (div-conforming) elements, it should be true.
%    bnd_sides: boundary sides on which the Dirichlet condition is imposed
%    bnd_func:  the condition to be imposed (g in the equations)
%    coeff:     function handle for the viscosity coefficient (mu in the equation)
%    Cpen:      a penalization term
%   
% OUTPUT:
%
%     vel:         assigned value to the normal degrees of freedom
%     normal_dofs: global numbering of the normal basis functions
%
% Copyright (C) 2014 Rafael Vazquez
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

function [u, dofs] = sp_drchlt_l2_proj_udotn (space, msh, bnd_sides, bnd_func)

  ndim = numel (msh.qn);

% The normal boundary condition is imposed strongly
  M = spalloc (space.ndof, space.ndof, 3*space.ndof);
  rhs2 = zeros (space.ndof, 1);

  normal_dofs = [];
  for iside = bnd_sides  
    msh_side = msh_eval_boundary_side (msh, iside);
    sp_side = sp_eval_boundary_side (space, msh_side);

    ind = floor ((iside+1)/2); % ind = [1 1 2 2];
    normal_dofs = union (normal_dofs, space.boundary(iside).comp_dofs{ind});

    for idim = 1:ndim
      x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    if (ndim == 2)
      g_func = @(x, y) bnd_func (x, y, iside);
    elseif (ndim == 3)
      g_func = @(x, y, z) bnd_func (x, y, z, iside);
    end
    g = g_func(x{:});
    
    M_loc = op_udotn_vdotn (sp_side, sp_side, msh_side, ones(size(x{1})));
    rhs_loc = op_fdotn_vdotn (sp_side, msh_side, g);
    
    M(sp_side.dofs, sp_side.dofs) = M(sp_side.dofs, sp_side.dofs) + M_loc;
    rhs2(sp_side.dofs) = rhs2(sp_side.dofs) + rhs_loc;
  end
  
  u = M(normal_dofs, normal_dofs) \ rhs2(normal_dofs);
  dofs = normal_dofs;
  
end
