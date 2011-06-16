% MSH_EVALUATE_COL: evaluate the parameterization in one column of the mesh.
%
%     msh = msh_evaluate_col (msh, geo)
%
% INPUTS:
%
%     msh:  mesh class (see msh_3d)
%     geo:  structure representing the geometrical mapping
%     'option', value: additional optional parameters, currently available options are:
%
% OUTPUT:
%
%     msh: structure containing the quadrature rule in one column of the physical domain, which contains the following fields
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     colnum        (scalar)                number of the column
%     nel           (scalar)                number of elements in the column
%     nelu          (scalar)                number of elements in the first parametric direction for the entire mesh
%     nelv          (scalar)                number of elements in the second parametric direction for the entire mesh
%     nelw          (scalar)                number of elements in the third parametric direction for the entire mesh
%     nqn           (scalar)                number of quadrature nodes per element
%     nqnu          (scalar)                number of quadrature nodes per element in the first parametric direction
%     nqnv          (scalar)                number of quadrature nodes per element in the second parametric direction
%     nqnw          (scalar)                number of quadrature nodes per element in the third parametric direction
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

  indu = colnum * ones(msh.nelv, msh.nelw);
  indv = repmat ((1:msh.nelv)', 1, msh.nelw);
  indw = repmat ((1:msh.nelw), msh.nelv, 1);
  elem_list = sub2ind ([msh.nelu, msh.nelv, msh.nelw], indu, indv, indw);
  msh_col.elem_list = elem_list(:);

  msh_col.nelu = msh.nelu;
  msh_col.nelv = msh.nelv;
  msh_col.nelw = msh.nelw;
  msh_col.nel  = msh.nelcol;

  msh_col.nqnu = msh.nqnu;
  msh_col.nqnv = msh.nqnv;
  msh_col.nqnw = msh.nqnw;
  msh_col.nqn  = msh.nqn;

  qnu = msh.qn{1}(:,colnum);  qnv = msh.qn{2}; qnw = msh.qn{3};

  quad_nodes_u = reshape (qnu, msh.nqnu, 1, 1, 1, 1);
  quad_nodes_u = repmat  (quad_nodes_u, [1, msh.nqnv, msh.nqnw, msh.nelv, msh.nelw]);
  quad_nodes_u = reshape (quad_nodes_u, [], msh_col.nel);

  quad_nodes_v = reshape (qnv, 1, msh.nqnv, msh.nelv, 1, 1);
  quad_nodes_v = repmat  (quad_nodes_v, [msh.nqnu, 1, msh.nqnw, 1, msh.nelw]);
  quad_nodes_v = reshape (quad_nodes_v, [], msh_col.nel);

  quad_nodes_w = reshape (qnw, 1, 1, 1, msh.nqnw, msh.nelw);
  quad_nodes_w = repmat  (quad_nodes_w, [msh.nqnu, msh.nqnv, 1, msh.nelv, 1]);
  quad_nodes_w = reshape (quad_nodes_w, [], msh_col.nel);

  msh_col.quad_nodes(1, :, :) = quad_nodes_u;
  msh_col.quad_nodes(2, :, :) = quad_nodes_v;
  msh_col.quad_nodes(3, :, :) = quad_nodes_w;

  clear quad_nodes_u quad_nodes_v quad_nodes_w

  if (~isempty (msh.qw))
    qwu = msh.qw{1}(:,colnum);  qwv = msh.qw{2}; qww = msh.qw{3};
    quad_weights_u = reshape (qwu, msh.nqnu, 1, 1, 1, 1);
    quad_weights_u = repmat  (quad_weights_u, [1, msh.nqnv, msh.nqnw, msh.nelv, msh.nelw]);
    quad_weights_u = reshape (quad_weights_u, [], msh_col.nel);

    quad_weights_v = reshape (qwv, 1, msh.nqnv, msh.nelv, 1, 1);
    quad_weights_v = repmat  (quad_weights_v, [msh.nqnu, 1, msh.nqnw, 1, msh.nelw]);
    quad_weights_v = reshape (quad_weights_v, [], msh_col.nel);

    quad_weights_w = reshape (qww, 1, 1, msh.nqnw, 1, msh.nelw);
    quad_weights_w = repmat  (quad_weights_w, [msh.nqnu, msh.nqnv, 1, msh.nelv, 1]);
    quad_weights_w = reshape (quad_weights_w, [], msh_col.nel);

    msh_col.quad_weights = quad_weights_u .* quad_weights_v .* quad_weights_w;

    clear quad_weights_u quad_weights_v quad_weights_w
  end

  F = feval (msh.map, {qnu(:)', qnv(:)', qnw(:)'});
  F = reshape (F, [3, msh.nqnu, msh.nqnv, msh.nelv, msh.nqnw, msh.nelw]);
  F = permute (F, [1 2 3 5 4 6]);
  msh_col.geo_map = reshape (F, [3, msh.nqn, msh_col.nel]);
  
  jac = feval (msh.map_der, {qnu(:)', qnv(:)' qnw(:)'});
  jac = reshape (jac, [3, 3, msh.nqnu, msh.nqnv, msh.nelv, msh.nqnw, msh.nelw]);
  jac = permute (jac, [1 2 3 4 6 5 7]);
  msh_col.geo_map_jac = reshape (jac, 3, 3, msh.nqn, msh_col.nel);
  msh_col.jacdet = abs (geopdes_det__ (msh_col.geo_map_jac));
  msh_col.jacdet = reshape (msh_col.jacdet, [msh_col.nqn, msh_col.nel]);

end
