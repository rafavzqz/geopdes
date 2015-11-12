% MP_SP_DRCHLT_L2_PROJ_UDOTN: assign the normal degrees of freedom trough an L2 projection for a multipatch geometry. 
%  To be used with the 'RT' and 'NDL' spaces. The imposed condition reads   u \cdot n = h \cdot n
%
%   [vel, normal_dofs] = mp_sp_drchlt_l2_proj_udotn (space, msh, gnum, ornt, bnd_sides, bnd_func)
%
% INPUTS:
%     
%    space:     space object (see sp_vector_div_transform)
%    msh:       mesh object (see msh_cartesian)
%    gnum:      global numbering of the degrees of freedom (see mp_interface_hdiv)
%    ornt:      global orientation of the degrees of freedom (see mp_interface_hdiv)
%  boundaries: array of structures containing the information for the boundaries (see mp_geo_load)
%    bnd_sides: boundary sides on which the Dirichlet condition is imposed
%    bnd_func:  the condition to be imposed (h in the equation)
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

function [u, dofs] = mp_sp_drchlt_l2_proj_udotn (space, msh, gnum, ornt, boundaries, refs, bnd_func)

  ndof = max ([gnum{:}]);
  M = spalloc (ndof, ndof, 3*ndof);
  rhs2 = zeros (ndof, 1);

  normal_dofs = [];
  for iref = refs
    for bnd_side = 1:boundaries(iref).nsides
      iptc = boundaries(iref).patches(bnd_side);
      iside = boundaries(iref).faces(bnd_side);

      if (isa (space{iptc}, 'sp_vector_div_transform'))
        [msh_side, msh_side_from_interior] = msh_eval_boundary_side (msh{iptc}, iside);
        sp_side = sp_eval_boundary_side (space{iptc}, msh_side, msh_side_from_interior);
      else
        error ('The function only works with div-conforming spaces')
      end

      ind = ceil (iside/2); % ind = [1 1 2 2] in 2D; ind = [1 1 2 2 3 3] in 3D
      normal_dofs = union (normal_dofs, gnum{iptc}(sp_side.comp_dofs{ind}));

      for idim = 1:msh{iptc}.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end
      g = bnd_func (x{:}, iref);

      M_loc = op_udotn_vdotn (sp_side, sp_side, msh_side, ones(size(x{1})));
      rhs_loc = op_fdotn_vdotn (sp_side, msh_side, g);
    
      global_dofs = gnum{iptc}(sp_side.dofs);
      M(global_dofs, global_dofs) = M(global_dofs, global_dofs) + M_loc;
      rhs2(global_dofs) = rhs2(global_dofs) + rhs_loc;
    end
  end
  
  u = M(normal_dofs, normal_dofs) \ rhs2(normal_dofs);
  dofs = normal_dofs;
  
end
