% MSH_EVALUATE_COL: evaluate the parameterization in one column of the mesh.
%
%     msh_col = msh_evaluate_col (msh, colnum)
%
% INPUTS:
%
%    msh:    mesh object (see msh_3d)
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
%     nel_dir       (1 x 3 vector)          number of elements in each parametric direction for the entire mesh
%     nqn           (scalar)                number of quadrature nodes per element
%     nqn_dir       (1 x 3 vector)          number of quadrature nodes per element in each parametric direction for the entire mesh
%     quad_nodes    (3 x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%     geo_map       (3 x nqn x nel vector)  physical coordinates of the quadrature nodes
%     geo_map_jac   (3 x 3 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
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

function msh_col = msh_evaluate_col (msh, colnum, varargin)

  msh_col.colnum = colnum;

  indu = colnum * ones(msh.nel_dir(2), msh.nel_dir(3));
  indv = repmat ((1:msh.nel_dir(2))', 1, msh.nel_dir(3));
  indw = repmat ((1:msh.nel_dir(3)), msh.nel_dir(2), 1);
  elem_list = sub2ind ([msh.nel_dir(1), msh.nel_dir(2), msh.nel_dir(3)], indu, indv, indw);
  msh_col.elem_list = elem_list(:);

  msh_col.nel_dir = msh.nel_dir;
  msh_col.nel  = msh.nelcol;

  msh_col.nqn_dir = msh.nqn_dir;
  msh_col.nqn  = msh.nqn;

  qnu = msh.qn{1}(:,colnum);  qnv = msh.qn{2}; qnw = msh.qn{3};

  if (isempty (msh.quad_nodes))
    quad_nodes_u = reshape (qnu, msh.nqn_dir(1), 1, 1, 1, 1, 1);
    quad_nodes_u = repmat  (quad_nodes_u, [1, msh.nqn_dir(2), msh.nqn_dir(3), 1, msh.nel_dir(2), msh.nel_dir(3)]);
    quad_nodes_u = reshape (quad_nodes_u, [], msh_col.nel);

    quad_nodes_v = reshape (qnv, 1, 1, msh.nqn_dir(2), msh.nel_dir(2), 1, 1);
    quad_nodes_v = repmat  (quad_nodes_v, [msh.nqn_dir(1), 1, msh.nqn_dir(3), 1, 1, msh.nel_dir(3)]);
    quad_nodes_v = reshape (quad_nodes_v, [], msh_col.nel);
  
    quad_nodes_w = reshape (qnw, 1, 1, 1, msh.nqn_dir(3), 1, msh.nel_dir(3));
    quad_nodes_w = repmat  (quad_nodes_w, [msh.nqn_dir(1), msh.nqn_dir(2), 1, 1, msh.nel_dir(2), 1]);
    quad_nodes_w = reshape (quad_nodes_w, [], msh_col.nel);

    msh_col.quad_nodes(1, :, :) = quad_nodes_u;
    msh_col.quad_nodes(2, :, :) = quad_nodes_v;
    msh_col.quad_nodes(3, :, :) = quad_nodes_w;

    clear quad_nodes_u quad_nodes_v quad_nodes_w
  else
    msh_col.quad_nodes = msh.quad_nodes(:,:,msh_col.elem_list);
  end


  if (~isempty (msh.qw))
    if (isempty (msh.quad_weights))
      qwu = msh.qw{1}(:,colnum);  qwv = msh.qw{2}; qww = msh.qw{3};
      quad_weights_u = reshape (qwu, msh.nqn_dir(1), 1, 1, 1, 1, 1);
      quad_weights_u = repmat  (quad_weights_u, [1, msh.nqn_dir(2), msh.nqn_dir(3), 1, msh.nel_dir(2), msh.nel_dir(3)]);
      quad_weights_u = reshape (quad_weights_u, [], msh_col.nel);

      quad_weights_v = reshape (qwv, 1, 1, msh.nqn_dir(2), msh.nel_dir(2), 1, 1);
      quad_weights_v = repmat  (quad_weights_v, [msh.nqn_dir(1), 1, msh.nqn_dir(3), 1, 1, msh.nel_dir(3)]);
      quad_weights_v = reshape (quad_weights_v, [], msh_col.nel);

      quad_weights_w = reshape (qww, 1, 1, msh.nqn_dir(3), 1, 1, msh.nel_dir(3));
      quad_weights_w = repmat  (quad_weights_w, [msh.nqn_dir(1), msh.nqn_dir(2), 1, 1, msh.nel_dir(2), 1]);
      quad_weights_w = reshape (quad_weights_w, [], msh_col.nel);

      msh_col.quad_weights = quad_weights_u .* quad_weights_v .* quad_weights_w;

      clear quad_weights_u quad_weights_v quad_weights_w
    else
      msh_col.quad_weights = msh.quad_weights(:,msh_col.elem_list);
    end
  end

  if (isempty (msh.geo_map))
    F = feval (msh.map, {qnu(:)', qnv(:)', qnw(:)'});
    F = reshape (F, [3, msh.nqn_dir(1), msh.nqn_dir(2), msh.nel_dir(2), msh.nqn_dir(3), msh.nel_dir(3)]);
    F = permute (F, [1 2 3 5 4 6]);
    msh_col.geo_map = reshape (F, [3, msh.nqn, msh_col.nel]);
  else
    msh_col.geo_map = msh.geo_map(:,:,msh_col.elem_list);
  end
  
  if (isempty (msh.geo_map_jac))
    jac = feval (msh.map_der, {qnu(:)', qnv(:)' qnw(:)'});
    jac = reshape (jac, [3, 3, msh.nqn_dir(1), msh.nqn_dir(2), msh.nel_dir(2), msh.nqn_dir(3), msh.nel_dir(3)]);
    jac = permute (jac, [1 2 3 4 6 5 7]);
    msh_col.geo_map_jac = reshape (jac, 3, 3, msh.nqn, msh_col.nel);
  else
    msh_col.geo_map_jac = msh.geo_map_jac(:,:,:,msh_col.elem_list);
  end

  if (isempty (msh.jacdet))
    msh_col.jacdet = abs (geopdes_det__ (msh_col.geo_map_jac));
    msh_col.jacdet = reshape (msh_col.jacdet, [msh_col.nqn, msh_col.nel]);
  else
    msh_col.jacdet = msh.jacdet(:,msh_col.elem_list);
  end

end
