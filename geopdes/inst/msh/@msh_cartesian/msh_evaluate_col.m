% MSH_EVALUATE_COL: evaluate the parameterization in one column of the mesh.
%
%     msh_col = msh_evaluate_col (msh, colnum)
%
% INPUTS:
%
%    msh:    mesh object (see msh_cartesian)
%    colnum: number of the "column", i.e., the element in the first parametric direction.
%
% OUTPUT:
%
%     msh_col: structure containing the quadrature rule in one column of the physical domain, which contains the following fields
%
%     FIELD_NAME    (SIZE)                    DESCRIPTION
%     ndim          (scalar)                  dimension of the parametric space
%     rdim          (scalar)                  dimension of the physical space
%     colnum        (scalar)                  number of the column
%     nel           (scalar)                  number of elements in the column
%     elem_list     (nel vector)              indices of the elements in the column
%     nel_dir       (1 x ndim vector)         number of elements in each parametric direction for the entire mesh
%     nqn           (scalar)                  number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)         number of quadrature nodes per element in each parametric direction for the entire mesh
%     quad_nodes    (ndim x nqn x nel vector) coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)        weights associated to the quadrature nodes
%     geo_map       (rdim x nqn x nel vector) physical coordinates of the quadrature nodes
%     geo_map_jac   (rdim x ndim x nqn x nel) Jacobian matrix of the map evaluated at the quadrature nodes
%     jacdet        (nqn x nel)               element of length, area, volume (if rdim = ndim, determinant of the Jacobian)
%     geo_map_der2  (rdim x ndim x ndim x nqn x nel) Hessian matrix of the map evaluated at the quadrature nodes
%     element_size  (1 x nel)                 characteristic size of the elements
%     normal        (rdim x ndim x nqn x nel) for 3D surfaces, the exterior normal to the surface
%
%  For more details, see the documentation
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Carlo de Falco, Sara Frizziero
% Copyright (C) 2015 Rafael Vazquez
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

  msh_col.ndim = msh.ndim;
  msh_col.rdim = msh.rdim;

  msh_col.colnum = colnum;

  if (msh.ndim == 1)
    msh_col.elem_list = colnum;
  elseif (msh.ndim == 2)
    msh_col.elem_list = colnum + msh.nel_dir(1)*(0:msh.nel_dir(2)-1);
  elseif (msh.ndim == 3)
    indu = colnum * ones(msh.nel_dir(2), msh.nel_dir(3));
    indv = repmat ((1:msh.nel_dir(2))', 1, msh.nel_dir(3));
    indw = repmat ((1:msh.nel_dir(3)), msh.nel_dir(2), 1);
    elem_list = sub2ind ([msh.nel_dir(1), msh.nel_dir(2), msh.nel_dir(3)], indu, indv, indw);
    msh_col.elem_list = elem_list(:);
  end

  msh_col.nel_dir = msh.nel_dir;
  msh_col.nel_dir(1) = 1;
  msh_col.nel = msh.nel / msh.nel_dir(1);
%   msh_col.nel  = msh.nelcol;

  msh_col.nqn_dir = msh.nqn_dir;
  msh_col.nqn  = msh.nqn;

  msh_col.qn = msh.qn; 
  msh_col.qn{1} = msh.qn{1}(:,colnum);
  if (~isempty (msh.qw))
    msh_col.qw = msh.qw; 
    msh_col.qw{1} = msh.qw{1}(:,colnum);
  end

  for idim = 1:msh.ndim
    qsize = ones (1, msh.ndim*2);
    qsize(2*idim-1:2*idim) = [msh_col.nqn_dir(idim), msh_col.nel_dir(idim)];
    qrep = [msh_col.nqn_dir(1:msh.ndim), msh_col.nel_dir(1:msh.ndim)];
    qrep([idim, msh.ndim+idim]) = 1;
    quad_nodes = reshape (msh_col.qn{idim}, qsize);
    quad_nodes = repmat (quad_nodes, qrep);
    msh_col.quad_nodes(idim,:,:) = reshape (quad_nodes, msh_col.nqn, msh_col.nel);
    clear qsize qrep quad_nodes
  end


  if (~isempty (msh.qw))
    qw = 1;
    for idim = 1:msh.ndim
      qw = kron (msh_col.qw{idim}, qw);
    end
    msh_col.quad_weights = qw;
  end
  
% Auxiliary vector sizes, to use with reshape and permute
  reorder = @(x) x(:)';
  psize = reorder ([msh_col.nqn_dir; msh_col.nel_dir]);
  vorder = [1:2:msh.ndim*2, 2:2:msh.ndim*2]; % [1 3 5 2 4 6], for ndim = 3
  F = feval (msh.map, cellfun (reorder, msh_col.qn, 'UniformOutput', false));
  F = reshape (F, [msh.rdim, psize]);
  F = permute (F, [1, vorder+1]);
  msh_col.geo_map = reshape (F, [msh.rdim, msh.nqn, msh_col.nel]);
  
  jac = feval (msh.map_der, cellfun (reorder, msh_col.qn, 'UniformOutput', false));
  jac = reshape (jac, [msh.rdim, msh.ndim, psize]);
  jac = permute (jac, [1 2 vorder+2]);
  msh_col.geo_map_jac = reshape (jac, msh.rdim, msh.ndim, msh.nqn, msh_col.nel);

  msh_col.jacdet = abs (geopdes_det__ (msh_col.geo_map_jac));
  msh_col.jacdet = reshape (msh_col.jacdet, [msh_col.nqn, msh_col.nel]);

  if (msh.der2)
    msh_col.geo_map_der2 = feval (msh.map_der2, cellfun (reorder, msh_col.qn, 'UniformOutput', false));
    msh_col.geo_map_der2 = reshape (msh_col.geo_map_der2, [msh.rdim, msh.ndim, msh.ndim, psize]);
    msh_col.geo_map_der2 = permute (msh_col.geo_map_der2, [1 2 3 vorder+3]);
    msh_col.geo_map_der2 = reshape (msh_col.geo_map_der2, [msh.rdim, msh.ndim, msh.ndim, msh.nqn, msh_col.nel]);
  end
  
  if (isfield(msh_col, 'quad_weights') && isfield(msh_col, 'jacdet'))
    msh_col.element_size = (sum (msh_col.quad_weights .* ...
                             abs (msh_col.jacdet), 1)).^(1/msh.ndim);
  end

  if (msh.ndim == 2 && msh.rdim == 3)
    normal = reshape (geopdes_cross__ (msh_col.geo_map_jac(:,1,:,:), ...
                              msh_col.geo_map_jac(:,2,:,:)), msh_col.rdim, msh_col.nqn, msh_col.nel);
    norms = reshape (geopdes_norm__ (normal), [1, msh_col.nqn, msh_col.nel]);
    msh_col.normal = bsxfun (@rdivide, normal, norms);
  end

end
