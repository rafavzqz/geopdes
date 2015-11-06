% MSH_EVALUATE_COL: evaluate the parameterization in one column of the mesh.
%
%     msh_col = msh_evaluate_col (msh, colnum)
%
% INPUTS:
%
%    msh:    mesh object (see msh_2d)
%    colnum: number of the "column", i.e., the element in the first parametric direction.
%
% OUTPUT:
%
%     msh_col: structure containing the quadrature rule in one column of the physical domain, which contains the following fields
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     colnum        (scalar)                number of the column
%     nel           (scalar)                number of elements in the column
%     elem_list     (nel vector)            indices of the elements in the column
%     nel_dir       (1 x 2 vector)          number of elements in each parametric direction for the entire mesh
%     nqn           (scalar)                number of quadrature nodes per element
%     nqn_dir       (1 x 2 vector)          number of quadrature nodes per element in each parametric direction
%     quad_nodes    (2 x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%     geo_map       (2 x nqn x nel vector)  physical coordinates of the quadrature nodes
%     geo_map_jac   (2 x 2 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%     geo_map_der2  (2 x 2 x 2 x nqn x nel) Second order derivatives of the map evaluated at the quadrature nodes (optional)
%     jacdet        (nqn x nel)             determinant of the Jacobian evaluated in the quadrature points
%  For more details, see the documentation
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

function msh_col = msh_evaluate_col (msh, colnum)

  msh_col.colnum = colnum;
  msh_col.elem_list = colnum + msh.nel_dir(1)*(0:msh.nel_dir(2)-1);

  msh_col.nel_dir = msh.nel_dir;
  msh_col.nel  = msh.nelcol;

  msh_col.nqn     = msh.nqn;
  msh_col.nqn_dir = msh.nqn_dir;

  qnu = msh.qn{1}(:,colnum);  qnv = msh.qn{2};

  if (isempty (msh.quad_nodes))
    quad_nodes_u = reshape (qnu, msh.nqn_dir(1), 1, 1);
    quad_nodes_u = repmat  (quad_nodes_u, [1, msh.nqn_dir(2), msh.nel_dir(2)]);
    quad_nodes_u = reshape (quad_nodes_u, [], msh.nel_dir(2));

    quad_nodes_v = reshape (qnv, 1, msh.nqn_dir(2), msh.nel_dir(2));
    quad_nodes_v = repmat  (quad_nodes_v, [msh.nqn_dir(1), 1, 1]);
    quad_nodes_v = reshape (quad_nodes_v, [], msh.nel_dir(2));

    msh_col.quad_nodes(1, :, :) = quad_nodes_u;
    msh_col.quad_nodes(2, :, :) = quad_nodes_v;

    clear quad_nodes_u quad_nodes_v
  else
    msh_col.quad_nodes = msh.quad_nodes(:,:,msh_col.elem_list);
  end


  if (~isempty (msh.qw))
    if (isempty (msh.quad_weights))
      qwu = msh.qw{1}(:,colnum);  qwv = msh.qw{2};
      quad_weights_u = reshape (qwu, msh.nqn_dir(1), 1, 1);
      quad_weights_u = repmat  (quad_weights_u, [1, msh.nqn_dir(2), msh.nel_dir(2)]);
      quad_weights_u = reshape (quad_weights_u, [], msh.nel_dir(2));

      quad_weights_v = reshape (qwv, 1, msh.nqn_dir(2), msh.nel_dir(2));
      quad_weights_v = repmat  (quad_weights_v, [msh.nqn_dir(1), 1, 1]);
      quad_weights_v = reshape (quad_weights_v, [], msh.nel_dir(2));

      msh_col.quad_weights = quad_weights_u .* quad_weights_v;

      clear quad_weights_u quad_weights_v
    else
      msh_col.quad_weights = msh.quad_weights(:,msh_col.elem_list);
    end
  end

  if (isempty (msh.geo_map))
    F = feval (msh.map, {qnu(:)', qnv(:)'});
    msh_col.geo_map = reshape (F, [2, msh.nqn, msh.nel_dir(2)]);
  else
    msh_col.geo_map = msh.geo_map(:,:,msh_col.elem_list);
  end
  
  if (isempty (msh.geo_map_jac))
    jac = feval (msh.map_der, {qnu(:)', qnv(:)'});
    msh_col.geo_map_jac = reshape (jac, 2, 2, msh.nqn, msh.nel_dir(2));
  else
    msh_col.geo_map_jac = msh.geo_map_jac(:,:,:,msh_col.elem_list);
  end
   
  if (isempty (msh.jacdet))
    msh_col.jacdet = abs (geopdes_det__ (msh_col.geo_map_jac));
    msh_col.jacdet = reshape (msh_col.jacdet, [msh.nqn, msh.nel_dir(2)]);
  else
    msh_col.jacdet = msh.jacdet(:,msh_col.elem_list);
  end

  if (msh.der2)
    if (isempty (msh.geo_map_der2))
      msh_col.geo_map_der2 = reshape (feval (msh.map_der2, {qnu(:)', qnv(:)'}), ...
                                         2, 2, 2, msh_col.nqn, msh_col.nel);
    else
      msh_col.geo_map_der2 = msh.geo_map_der2(:,:,:,:,msh_col.elem_list);
    end
  end

end
