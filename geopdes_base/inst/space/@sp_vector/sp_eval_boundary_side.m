% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary.
%
%     sp_side = sp_eval_boundary_side (sp, msh_side)
%
% INPUTS:
%
%     sp:       space object (see sp_vector)
%     msh_side: mesh structure containing the information of the quadrature
%               rule on the boundary edge (see msh_cartesian/msh_eval_boundary_side)
%
% OUTPUT:
%
%     sp_side: structure that contains the following fields
%              (see the article for a detailed description)
%
%     FIELD_NAME      (SIZE)                                          DESCRIPTION
%     ncomp           (scalar)                                        number of components of the functions of the space (actually, 3)
%     nsh_max         (scalar)                                        maximum number of shape functions per element
%     nsh             (1 x msh_side.nel vector)                       actual number of shape functions per each element
%     ndof            (scalar)                                        total number of degrees of freedom
%     connectivity    (nsh_max x msh_side.nel vector)                 indices of basis functions that do not vanish in each element
%     shape_functions (ncomp x msh_side.nqn x nsh_max x msh_side.nel) basis functions evaluated at each quadrature node in each element
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
  sp_bnd = sp.boundary(iside);

  
  sp_bnd.nsh = zeros (1, msh_side.nel);
  sp_bnd.connectivity = [];
  sp_bnd.shape_functions = zeros (sp.ncomp, msh_side.nqn, sp_bnd.nsh_max, msh_side.nel);

  aux = 0; aux2 = 0;
  for icomp = 1:sp.ncomp
    sp_bnd_scalar = sp_eval_boundary_side (sp.scalar_spaces{icomp}, msh_side);
    sp_bnd.nsh = sp_bnd.nsh + sp_bnd_scalar.nsh;
    sp_bnd.connectivity = [sp_bnd.connectivity; aux+sp_bnd_scalar.connectivity];
    sp_bnd.shape_functions(icomp,:,aux2+(1:sp_bnd_scalar.nsh_max),:) = sp_bnd_scalar.shape_functions;
    aux  = aux + sp_bnd_scalar.ndof;
    aux2 = aux2 + sp_bnd_scalar.nsh_max;
  end

end