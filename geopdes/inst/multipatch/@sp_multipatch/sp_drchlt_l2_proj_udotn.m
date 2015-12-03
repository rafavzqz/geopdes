% SP_DRCHLT_L2_PROJ_UDOTN: assign the normal degrees of freedom through an L2 projection for a multipatch geometry. 
%  To be used with the 'RT' and 'NDL' spaces (div-preserving). The imposed condition reads   u \cdot n = h \cdot n
%
%   [vel, normal_dofs] = sp_drchlt_l2_proj_udotn (space, msh, bnd_sides, bnd_func)
%
% INPUTS:
%     
%    space:      multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch). 
%    msh:        multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%    bnd_sides:  boundary sides on which the Dirichlet condition is imposed
%    bnd_func:   the condition to be imposed (h in the equation)
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

function [u, normal_dofs] = sp_drchlt_l2_proj_udotn (space, msh, refs, bnd_func, varargin)

  M = spalloc (space.boundary.ndof, space.boundary.ndof, 3*space.boundary.ndof);
  rhs = zeros (space.boundary.ndof, 1);
  
  boundaries = msh.boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  bnd_dofs = [];
  for iref = refs
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    href = @(varargin) bnd_func(varargin{:}, iref);
    f_one = @(varargin) ones (size(varargin{1}));
      
    boundary_gnum = space.boundary.gnum;
    
    M = M + op_u_v_mp (space.boundary, space.boundary, msh.boundary, f_one, iref_patch_list);
    
% For this part we need the normal, which is computed with eval_boundary_side    
    for iptc = iref_patch_list
      patch_number = msh.boundary.patch_numbers(iptc);
      side_number = msh.boundary.side_numbers(iptc);
      msh_side = msh_eval_boundary_side (msh.msh_patch{patch_number}, side_number);
      sp_side = sp_eval_boundary_side (space.sp_patch{patch_number}, msh_side);
      
      for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
      end
      rhs(boundary_gnum{iptc}) = rhs(boundary_gnum{iptc}) + op_fdotn_v (sp_side, msh_side, href(x{:}));
    end
    
    boundary_gnum = space.boundary.gnum;
    bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}]);
  end
    
  u = M(bnd_dofs,bnd_dofs) \ rhs(bnd_dofs);
  normal_dofs = space.boundary.dofs(bnd_dofs);
  
end
