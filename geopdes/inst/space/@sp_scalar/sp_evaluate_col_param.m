% SP_EVALUATE_COL_PARAM: compute the basis functions, in the parametric domain, in one column of the mesh.
%
%     sp = sp_evaluate_col_param (space, msh_col, 'option1', value1, ...)
%
% INPUTS:
%
%    space:   object defining the space of discrete functions (see sp_scalar)
%    msh_col: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_cartesian/msh_evaluate_col)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                 DESCRIPTION
%    ncomp           (scalar)                               number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                               total number of degrees of freedom
%    ndof_dir        (1 x ndim vector)                      degrees of freedom along each direction
%    nsh_max         (scalar)                               maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%         (ndim x msh_col.nqn x nsh_max x msh_col.nel)      basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%         (ndim x ndim x msh_col.nqn x nsh_max x msh_col.nel)  basis function hessians evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
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

function sp = sp_evaluate_col_param (space, msh, varargin)

value = true;
gradient = false;
hessian = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col_param: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    else
      error ('sp_evaluate_col_param: unknown option %s', varargin {ii});
    end
  end
end

sp_univ = space.sp_univ;

elem_list{1} = msh.colnum;
for idim = 2:msh.ndim
  elem_list{idim} = 1:msh.nel_dir(idim);
end

for idim = 1:msh.ndim
  nsh_dim{idim} = sp_univ(idim).nsh(elem_list{idim});
end

[nsh_grid{1:msh.ndim}] = ndgrid (nsh_dim{:});
nsh = 1;
for idim = 1:msh.ndim
  nsh = nsh .* nsh_grid{idim};
end
nsh = nsh(:)';

for idim = 1:msh.ndim
  csize = ones (1, 2*msh.ndim);
  csize([idim, msh.ndim+idim]) = [sp_univ(idim).nsh_max, msh.nel_dir(idim)];
  crep = [sp_univ.nsh_max, msh.nel_dir];
  crep([idim, msh.ndim+idim]) = 1;

  conn{idim} = reshape (sp_univ(idim).connectivity(:,elem_list{idim}), csize);
  conn{idim} = repmat (conn{idim}, crep);
  conn{idim} = reshape (conn{idim}, [], msh.nel);
end

connectivity = zeros (space.nsh_max, msh.nel);
indices = ones (size (conn{1}));
for idim = 1:msh.ndim
  indices = indices & conn{idim} ~= 0;
end
for idim = 1:msh.ndim
  conn{idim} = conn{idim}(indices);
end
connectivity(indices) = sub2ind ([space.ndof_dir, 1], conn{:}); % The extra 1 makes things work in any dimension
connectivity = reshape (connectivity, space.nsh_max, msh.nel);

clear conn csize crep indices

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', space.ndof,  ...
            'ndof_dir', space.ndof_dir, 'connectivity', connectivity, ...
            'ncomp', 1, 'degree', space.degree);

shp = cell(1,msh.ndim); shg = cell(1,msh.ndim); shh = cell(1,msh.ndim);
for idim = 1:msh.ndim
  ssize = ones (1, 3*msh.ndim);
  ssize([idim, msh.ndim+idim, 2*msh.ndim+idim]) = [msh.nqn_dir(idim), sp_univ(idim).nsh_max, msh.nel_dir(idim)];
  srep = [msh.nqn_dir, sp_univ.nsh_max, msh.nel_dir];
  srep([idim, msh.ndim+idim, 2*msh.ndim+idim]) = 1;
  shp{idim} = reshape (sp_univ(idim).shape_functions(:,:,elem_list{idim}), ssize);
  shp{idim} = repmat (shp{idim}, srep);
  shp{idim} = reshape (shp{idim}, msh.nqn, space.nsh_max, msh.nel);  
  shg{idim} = reshape (sp_univ(idim).shape_function_gradients(:,:,elem_list{idim}), ssize);
  shg{idim} = repmat (shg{idim}, srep);
  shg{idim} = reshape (shg{idim}, msh.nqn, space.nsh_max, msh.nel);  
  shh{idim} = reshape (sp_univ(idim).shape_function_hessians(:,:,elem_list{idim}), ssize);
  shh{idim} = repmat (shh{idim}, srep);
  shh{idim} = reshape (shh{idim}, msh.nqn, space.nsh_max, msh.nel);  
end

if (value)
  sp.shape_functions = 1;
  for idim = 1:msh.ndim
    sp.shape_functions = sp.shape_functions .* shp{idim};
  end
end

if (gradient)
  for idim = 1:msh.ndim
    shape_fun_grad = shg{idim};
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_grad = shape_fun_grad .* shp{jdim};
    end
    sp.shape_function_gradients(idim,:,:,:) = shape_fun_grad;
  end
  sp.shape_function_gradients = reshape (sp.shape_function_gradients, ...
                                msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
end

if (hessian && isfield (msh, 'geo_map_der2'))
  for idim = 1:msh.ndim
    shape_fun_hess = shh{idim};
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_hess = shape_fun_hess .* shp{jdim};
    end
    sp.shape_function_hessians(idim,idim,:,:,:) = shape_fun_hess;
    
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_hess = shg{idim} .* shg{jdim};
      for kdim = setdiff (1:msh.ndim, [idim, jdim])
        shape_fun_hess = shape_fun_hess .* shp{kdim};
      end
      sp.shape_function_hessians(idim,jdim,:,:,:) = shape_fun_hess;
    end
  end
end

clear shp shg shh

if (value || gradient || hessian)
  if (strcmpi (space.space_type, 'NURBS'))
    sp = bsp_2_nrb__ (sp, msh, space.weights);
  end
end
    
end
