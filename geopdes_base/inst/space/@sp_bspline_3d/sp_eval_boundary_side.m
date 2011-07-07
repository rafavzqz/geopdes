% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary.
%
%     sp_side = sp_eval_boundary_side (sp, msh_side)
%
% INPUTS:
%
%     sp:       space object (see sp_bspline_3d)
%     msh_side: mesh structure containing the information of the quadrature
%               rule on the boundary edge (see msh_3d/msh_eval_boundary_side)
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
  sp_bnd = sp.boundary(iside);

  ind = setdiff (1:3, ceil(iside/2)); % ind = [2 3; 2 3; 1 3; 1 3; 1 2;	1 2] 

% Set the two parametric directions of the boundary
  switch (iside)
    case {1, 2}
      spu = sp.spv;
      spv = sp.spw;
    case {3, 4}
      spu = sp.spu;
      spv = sp.spw;
    case {5, 6}
      spu = sp.spu;
      spv = sp.spv;
  end

  nsh  = spu.nsh' * spv.nsh;
  sp_bnd.nsh  = nsh(:)';

  nel_dir = msh_side.nel_dir;
  nqn_dir = msh_side.nqn_dir;
  nel = msh_side.nel;
  nqn = msh_side.nqn;

  conn_u = reshape (spu.connectivity, spu.nsh_max, 1, nel_dir(1), 1);
  conn_u = repmat  (conn_u, [1, spv.nsh_max, 1, nel_dir(2)]);
  conn_u = reshape (conn_u, [], nel);

  conn_v = reshape (spv.connectivity, 1, spv.nsh_max, 1, nel_dir(2));
  conn_v = repmat  (conn_v, [spu.nsh_max, 1, nel_dir(1), 1]);
  conn_v = reshape (conn_v, [], nel);

  connectivity = zeros (sp_bnd.nsh_max, nel);
  indices = (conn_u ~= 0) & (conn_v ~= 0);
  connectivity(indices) = ...
      sub2ind ([spu.ndof, spv.ndof], conn_u(indices), conn_v(indices));
  sp_bnd.connectivity = reshape (connectivity, sp_bnd.nsh_max, nel);

  clear conn_u conn_v

  shp_u = reshape (spu.shape_functions, nqn_dir(1), 1, spu.nsh_max, 1, nel_dir(1), 1);
  shp_u = repmat  (shp_u, [1, nqn_dir(2), 1, spv.nsh_max, 1, nel_dir(2)]);
  shp_u = reshape (shp_u, nqn, sp_bnd.nsh_max, nel);

  shp_v = reshape (spv.shape_functions, 1, nqn_dir(2), 1, spv.nsh_max, 1, nel_dir(2));
  shp_v = repmat  (shp_v, [nqn_dir(1), 1, spu.nsh_max, 1, nel_dir(1), 1]);
  shp_v = reshape (shp_v, nqn, sp_bnd.nsh_max, nel);

  sp_bnd.shape_functions = shp_u .* shp_v ;
  clear shp_u shp_v

end
