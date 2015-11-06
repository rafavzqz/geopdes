% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary, for tangential basis functions.
%
%     sp_side = sp_eval_boundary_side (sp, msh_side)
%
% INPUTS:
%
%     sp:       space object (see sp_vector_3d_curl_transform)
%     msh_side: mesh structure containing the information of the quadrature
%               rule on the boundary edge (see msh_3d/msh_eval_boundary_side)
%
% OUTPUT:
%
%     sp_side: structure that contains the following fields
%              (see the article for a detailed description)
%
%     FIELD_NAME      (SIZE)                                       DESCRIPTION
%     ncomp           (scalar)                                     number of components of the functions of the space (actually, 2)
%     nsh_max         (scalar)                                     maximum number of shape functions per element
%     nsh             (1 x msh_side.nel vector)                    actual number of shape functions per each element
%     ndof            (scalar)                                     total number of degrees of freedom
%     connectivity    (nsh_max x msh_side.nel vector)              indices of basis functions that do not vanish in each element
%     shape_functions (2 x msh_side.nqn x nsh_max x msh_side.nel)  tangential basis functions evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
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

  iface = msh_side.side_number;
  sp_bnd = sp.boundary(iface);

  sp_bnd1 = sp_eval_boundary_side (sp.sp1, msh_side);
  sp_bnd2 = sp_eval_boundary_side (sp.sp2, msh_side);
  sp_bnd3 = sp_eval_boundary_side (sp.sp3, msh_side);

% Fields for the transformation
  [JinvT, jacdet] = geopdes_invT__ (msh_side.geo_map_jac);
  JinvT = reshape (JinvT, [3, 3, msh_side.nqn, msh_side.nel]);

  switch (iface)
   case {1, 2}
     sp_bnd.nsh = sp_bnd2.nsh + sp_bnd3.nsh;
     sp_bnd.connectivity  = [sp_bnd2.connectivity; sp_bnd3.connectivity + sp_bnd2.ndof];

     sp2_shape_funs = reshape (sp_bnd2.shape_functions, ...
          [msh_side.nqn, sp_bnd2.nsh_max, msh_side.nel]);
     sp3_shape_funs = reshape (sp_bnd3.shape_functions, ...
          [msh_side.nqn, sp_bnd3.nsh_max, msh_side.nel]);

     sp_bnd.shape_functions = zeros (3, msh_side.nqn, sp_bnd.nsh_max, msh_side.nel);
     sp_bnd.shape_functions(2,:,1:sp_bnd2.nsh_max,:) = sp2_shape_funs;
     sp_bnd.shape_functions(3,:,(sp_bnd2.nsh_max+1):sp_bnd.nsh_max,:) = sp3_shape_funs;

     sp_bnd.shape_functions = geopdes_prod__ (JinvT, sp_bnd.shape_functions);
   case {3, 4}
     sp_bnd.nsh = sp_bnd1.nsh + sp_bnd3.nsh;
     sp_bnd.connectivity  = [sp_bnd1.connectivity; sp_bnd3.connectivity + sp_bnd1.ndof];

     sp1_shape_funs = reshape (sp_bnd1.shape_functions, ...
          [msh_side.nqn, sp_bnd1.nsh_max, msh_side.nel]);
     sp3_shape_funs = reshape (sp_bnd3.shape_functions, ...
          [msh_side.nqn, sp_bnd3.nsh_max, msh_side.nel]);

     sp_bnd.shape_functions = zeros (3, msh_side.nqn, sp_bnd.nsh_max, msh_side.nel);
     sp_bnd.shape_functions(1,:,1:sp_bnd1.nsh_max,:) = sp1_shape_funs;
     sp_bnd.shape_functions(3,:,(sp_bnd1.nsh_max+1):sp_bnd.nsh_max,:) = sp3_shape_funs;

     sp_bnd.shape_functions = geopdes_prod__ (JinvT, sp_bnd.shape_functions);
   case {5, 6}
     sp_bnd.nsh = sp_bnd1.nsh + sp_bnd2.nsh;
     sp_bnd.connectivity  = [sp_bnd1.connectivity; sp_bnd2.connectivity + sp_bnd1.ndof];

     sp1_shape_funs = reshape (sp_bnd1.shape_functions, ...
          [msh_side.nqn, sp_bnd1.nsh_max, msh_side.nel]);
     sp2_shape_funs = reshape (sp_bnd2.shape_functions, ...
          [msh_side.nqn, sp_bnd2.nsh_max, msh_side.nel]);

     sp_bnd.shape_functions = zeros (3, msh_side.nqn, sp_bnd.nsh_max, msh_side.nel);
     sp_bnd.shape_functions(1,:,1:sp_bnd1.nsh_max,:) = sp1_shape_funs;
     sp_bnd.shape_functions(2,:,(sp_bnd1.nsh_max+1):sp_bnd.nsh_max,:) = sp2_shape_funs;

     sp_bnd.shape_functions = geopdes_prod__ (JinvT, sp_bnd.shape_functions);
  end

end
