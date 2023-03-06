% MSH_EVALUATE_ELEMENT_LIST: evaluate the parameterization in a given list of elements.
%
%     msh_elems = msh_evaluate_element_list (msh, elements)
%
% INPUTS:
%
%    msh:          mesh object (see msh_cartesian)
%    element_list: numbering of the elements where the evaluations are performed.
%
% OUTPUT:
%
%     msh_elems: structure containing the quadrature rule in the given elements of the physical domain, which contains the following fields
%
%     FIELD_NAME    (SIZE)                    DESCRIPTION
%     ndim          (scalar)                  dimension of the parametric space
%     rdim          (scalar)                  dimension of the physical space
%     nel           (scalar)                  number of elements in the list
%     elem_list     (1 x nel)                 numbering of the elements in the list
%     nel_dir       (1 x ndim vector)         number of elements in each parametric direction for the entire mesh
%     nqn           (scalar)                  number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)         number of quadrature nodes per element in each parametric direction for the entire mesh
%     der2          (scalar)                  wether second derivatives must be computed
%     der3          (scalar)                  wether third derivatives must be computed
%     der4          (scalar)                  wether fourth derivatives must be computed
%     quad_nodes    (ndim x nqn x nel vector) coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)        weights associated to the quadrature nodes
%     geo_map       (rdim x nqn x nel vector) physical coordinates of the quadrature nodes
%     geo_map_jac   (rdim x ndim x nqn x nel) Jacobian matrix of the map evaluated at the quadrature nodes
%     jacdet        (nqn x nel)               element of length, area, volume (if rdim = ndim, determinant of the Jacobian)
%     geo_map_der2  (rdim x ndim x ndim x nqn x nel])
%                                             Hessian tensor of the map evaluated at the quadrature nodes
%     geo_map_der3  (rdim x ndim x ndim x ndim x nqn x nel])
%                                             Third derivatives tensor of the map evaluated at the quadrature nodes
%     geo_map_der4  (rdim x ndim x ndim x ndim x ndim x nqn x nel])
%                                             Fourth derivatives tensor of the map evaluated at the quadrature nodes
%     normal        (rdim x ndim x nqn x nel]) for 3D surfaces, the exterior normal to the surface
%
%  For more details, see the documentation
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2015 Rafael Vazquez
% Copyright (C) 2023 Pablo Antolin
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

function msh_col = msh_evaluate_element_list (msh, elem_list, varargin)

  elem_list = elem_list(:)';

  msh_col.ndim = msh.ndim;
  msh_col.rdim = msh.rdim;

  msh_col.nel  = numel (elem_list);
  msh_col.elem_list = elem_list;
  msh_col.nel_dir = msh.nel_dir;

  msh_col.nqn_dir = msh.nqn_dir;
  msh_col.nqn  = msh.nqn;

  msh_col.der2 = msh.der2;
  msh_col.der3 = msh.der3;
  msh_col.der4 = msh.der4;

  indices = cell (msh.ndim, 1);
  [indices{:}] = ind2sub (msh.nel_dir, elem_list);
  indices = cell2mat (indices);

  if (isempty (elem_list))
    msh_col.quad_weights = [];
    msh_col.geo_map = [];
    msh_col.geo_map_jac = [];
    msh_col.geo_map_der2 = [];
    msh_col.geo_map_der3 = [];
    msh_col.geo_map_der4 = [];
    msh_col.jacdet = [];
    msh_col.element_size = [];
    return
  end

  qn_elems = arrayfun(@(ii) {msh.qn{ii}(:,indices(ii,:))}, 1:msh.ndim);
  qqn = cell (1,msh_col.nel);
  for iel = 1:numel(elem_list)
    for idim = 1:msh.ndim
      qqn{iel}{idim} = qn_elems{idim}(:,iel)';
    end
  end

  for iel = 1:numel(elem_list)
    xx = cell (msh.ndim, 1);
    [xx{:}] = ndgrid (qqn{iel}{:});
    for idim = 1:msh.ndim
      quad_nodes(idim,:,iel) = xx{idim}(:)';
    end

    if (~isempty (msh.qw))
      qw = 1;
      for idim = 1:msh.ndim
        qw = kron (msh.qw{idim}(:,indices(idim,iel)), qw);
      end
      msh_col.quad_weights(:,iel) = qw;
    end
  end

  msh_col.geo_map = zeros (msh.rdim, msh.nqn, msh_col.nel);
  msh_col.geo_map_jac  = zeros (msh.rdim, msh.ndim, msh.nqn, msh_col.nel);
  if (msh.der2)
    msh_col.geo_map_der2 = zeros (msh.rdim, msh.ndim, msh.ndim, msh.nqn, msh_col.nel);
    if (msh.der3)
      msh_col.geo_map_der3 = zeros (msh.rdim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, msh_col.nel);
      if (msh.der4)
        msh_col.geo_map_der4 = zeros (msh.rdim, msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, msh_col.nel);
        [F, jac, hess, der3, der4] = cellfun (@(x) feval (msh.map_der4, x), qqn, 'UniformOutput', false);
      else
        [F, jac, hess, der3] = cellfun (@(x) feval (msh.map_der3, x), qqn, 'UniformOutput', false);
      end
    else
      [F, jac, hess] = cellfun (@(x) feval (msh.map_der2, x), qqn, 'UniformOutput', false);
    end
  else
    [F, jac] = cellfun (@(x) feval (msh.map_der, x), qqn, 'UniformOutput', false);
  end

  for iel = 1:numel(elem_list)
    msh_col.geo_map(:,:,iel) = reshape (F{iel}, msh.rdim, msh.nqn);
    msh_col.geo_map_jac(:,:,:,iel) = reshape (jac{iel}, msh.rdim, msh.ndim, msh.nqn);
  end

  if (msh.der2)
    for iel = 1:numel(elem_list)
      msh_col.geo_map_der2(:,:,:,:,iel) = reshape (hess{iel}, msh.rdim, msh.ndim, msh.ndim, msh.nqn);
    end
  end
  if (msh.der3)
    for iel = 1:numel(elem_list)
      msh_col.geo_map_der3(:,:,:,:,:,iel) = reshape (der3{iel}, msh.rdim, msh.ndim, msh.ndim, msh.ndim, msh.nqn);
    end
  end
  if (msh.der4)
    for iel = 1:numel(elem_list)
      msh_col.geo_map_der4(:,:,:,:,:,:,iel) = reshape (der4{iel}, msh.rdim, msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn);
    end
  end


  msh_col.jacdet = abs (geopdes_det__ (msh_col.geo_map_jac));
  msh_col.jacdet = reshape (msh_col.jacdet, [msh_col.nqn, msh_col.nel]);

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
