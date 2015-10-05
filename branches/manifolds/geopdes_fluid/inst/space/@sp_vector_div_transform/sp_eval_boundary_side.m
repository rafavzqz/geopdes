% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary.
%
%     sp_side = sp_eval_boundary_side (sp, msh_side)
%
% INPUTS:
%
%     sp:       space object (see sp_vector_div_transform)
%     msh_side: mesh structure containing the information of the quadrature
%               rule on the boundary edge (see msh_cartesian/msh_eval_boundary_side)
%
% OUTPUT:
%
%     sp_side: structure that contains the following fields
%              (see the article for a detailed description)
%
%     FIELD_NAME      (SIZE)                                           DESCRIPTION
%     ncomp           (scalar)                                         number of components of the functions of the space (actually, 3)
%     nsh_max         (scalar)                                         maximum number of shape functions per element
%     nsh             (1 x msh_side.nel vector)                        actual number of shape functions per each element
%     ndof            (scalar)                                         total number of degrees of freedom
%     connectivity    (nsh_max x msh_side.nel vector)                  indices of basis functions that do not vanish in each element
%     shape_functions (ncomp x msh_side.nqn x nsh_max x msh_side.nel)  basis functions evaluated at each quadrature node in each element
%
% Copyright (C) 2013 Elena Bulgarello, Carlo de Falco, Sara Frizziero
% Copyright (C) 2015 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function sp_bnd = sp_eval_boundary_side (sp, msh_side, msh_aux)

  iside = msh_side.side_number;
  sp_bnd = sp.boundary(iside);

  for icomp = 1:sp.ncomp_param
    sp_bnd_scalar{icomp} = sp_eval_boundary_side (sp.scalar_spaces{icomp}, msh_side);
  end

  nsh  = zeros (1, msh_side.nel);
  connectivity = [];
  aux = 0;
  for icomp = 1:sp.ncomp_param
    ndof_dir(icomp,:) = sp_bnd_scalar{icomp}.ndof_dir;
    nsh = nsh + sp_bnd_scalar{icomp}.nsh(:)';
  
    connectivity = [connectivity; sp_bnd_scalar{icomp}.connectivity+aux];
    aux = aux + sp_bnd_scalar{icomp}.ndof;
  end
  
  sp_bnd.nsh = nsh;
  sp_bnd.connectivity = connectivity;

  aux = 0;
  shape_funs = zeros (sp.ncomp_param, msh_side.nqn, sp_bnd.nsh_max, msh_side.nel);
  for icomp = 1:sp.ncomp_param
    indices = aux+(1:sp_bnd_scalar{icomp}.nsh_max);
    shape_funs(icomp,:,indices,:) = sp_bnd_scalar{icomp}.shape_functions;
    aux = indices(end);
  end
  
  shape_funs = geopdes_prod__ (msh_aux.geo_map_jac, shape_funs);
  jacdet = reshape (geopdes_det__ (msh_aux.geo_map_jac), 1, msh_side.nqn, 1, msh_side.nel);
  sp_bnd.shape_functions = bsxfun (@rdivide, shape_funs, jacdet);

end
