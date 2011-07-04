% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary.
%
%     sp_side = sp_eval_boundary_side (sp, msh_side)
%
% INPUTS:
%
%     sp:       space class (see sp_nurbs_2d)
%     msh_side: mesh structure containing the information of the quadrature
%               rule on the boundary edge (see msh_2d/msh_eval_boundary_side)
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

function sp_bnd = sp_eval_boundary_side (sp, msh_side)

  iside = msh_side.side_number;
  ind = mod (floor ((iside+1)/2), 2) + 1;

  bnodes = reshape (squeeze (msh_side.quad_nodes(ind,:,:)), msh_side.nqn, []);
  sp_bnd = sp_bspline_1d_param (sp.knots{ind}, sp.degree(ind), bnodes, ...
                                'gradient', false);

  switch (iside)
   case 1
    weights = sp.weights(1,:);
   case 2
    weights = sp.weights(end,:);
   case 3
    weights = sp.weights(:,1);
   case 4
    weights = sp.weights(:,end);
  end

  sp_bnd = bsp_2_nrb_1d__ (sp_bnd, msh_side, weights);

  sp_bnd.dofs = sp.boundary(iside).dofs;

end
