% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary.
%
%     sp_side = sp_eval_boundary_side (sp, msh_side)
%
% INPUTS:
%
%     sp:       space object (see sp_nurbs)
%     msh_side: mesh structure containing the information of the quadrature
%               rule on the boundary edge (see msh_cartesian/msh_eval_boundary_side)
%
% OUTPUT:
%
%     sp_side: structure that contains the following fields
%              (see the article for a detailed description)
%
%     FIELD_NAME      (SIZE)                                   DESCRIPTION
%     ncomp           (scalar)                                 number of components of the functions of the space (actually, 1)
%     nsh_max         (scalar)                                 maximum number of shape functions per element
%     nsh             (1 x msh_side.nel vector)                actual number of shape functions per each element
%     ndof            (scalar)                                 total number of degrees of freedom
%     connectivity    (nsh_max x msh_side.nel vector)          indices of basis functions that do not vanish in each element
%     shape_functions (msh_side.nqn x nsh_max x msh_side.nel)  basis functions evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function sp_bnd = sp_eval_boundary_side (sp, msh_side)

  iside = msh_side.side_number;
  sp_bnd = sp_eval_boundary_side (sp.spline_space, msh_side);

%%    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
%%    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
  ind2 = ceil (iside/2);
%   ind = setdiff (1:msh_side.ndim+1, ind2);

  for idim = 1:msh_side.ndim+1
    dofs_dir{idim} = 1:sp.ndof_dir(idim);
  end
  
  if (rem (iside, 2) == 1)
    dofs_dir{ind2} = 1;
  elseif (rem (iside, 2) == 0)
    dofs_dir{ind2} = sp.ndof_dir(ind2);
  end
  
  weights = sp.weights(dofs_dir{:});
  sp_bnd = bsp_2_nrb__ (sp_bnd, msh_side, weights);

end
