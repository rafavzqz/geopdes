% SP_EVALUATE_COL: compute the basis functions in one column of the mesh.
%
%     sp = sp_evaluate_col (space, msh_col, 'option1', value1, ...)
%
% INPUTS:
%     
%    sp:      object defining the space of discrete functions (see sp_vector_div_transform)
%    msh_col: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_cartesian/msh_evaluate_col)
%    option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            divergence |      false      |  compute shape_function_divs
%            curl       |      false      |  compute shape_function_curls
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                        DESCRIPTION
%    ncomp           (scalar)                                      number of components of the functions of the space (actually, 3)
%    ndof            (scalar)                                      total number of degrees of freedom
%    ndof_dir        (ncomp x ndim matrix)                         for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                                      maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)                      actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)                indices of basis functions that do not vanish in each element
%    shape_functions (ncomp x msh_col.nqn x nsh_max x msh_col.nel) basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%           (ncomp x rdim x msh_col.nqn x nsh_max x msh_col.nel)   basis function gradients evaluated at each quadrature node in each element
%    shape_function_divs (msh_col.nqn x nsh_max x msh_col.nel)     basis function gradients evaluated at each quadrature node in each element
%    shape_function_curls (msh_col.nqn x nsh_max x msh_col.nel)    basis function gradients evaluated at each quadrature node in each element
%
% Copyright (C) 2013 Elena Bulgarello, Carlo de Falco, Sara Frizziero
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

function sp = sp_evaluate_col (space, msh, varargin)

value = true;
gradient = false;
divergence = false;
curl = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'curl'))
      curl = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'divergence'))
      divergence = varargin {ii+1};
    else
      error ('sp_evaluate_col: unknown option %s', varargin {ii});
    end
  end
end

first_der = gradient || divergence || curl;
for icomp = 1:space.ncomp_param
  sp_col_scalar(icomp) = sp_evaluate_col_param (space.scalar_spaces{icomp}, msh, 'value', true, 'gradient', first_der);
end

ndof_scalar = [sp_col_scalar.ndof];

ndof = sum (ndof_scalar);
nsh  = zeros (1, msh.nel);
connectivity = [];
aux = 0;
for icomp = 1:space.ncomp_param
  ndof_dir(icomp,:) = sp_col_scalar(icomp).ndof_dir;
  nsh = nsh + sp_col_scalar(icomp).nsh(:)';
  
  connectivity = [connectivity; sp_col_scalar(icomp).connectivity+aux];
  aux = aux + ndof_scalar(icomp);
end


sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
            'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
            'ncomp', space.ncomp, 'ncomp_param', space.ncomp_param);

% We start computing the shape functions and gradients in the reference domain
shape_functions = zeros (sp.ncomp_param, msh.nqn, sp.nsh_max, msh.nel);
for icomp = 1:space.ncomp_param
  indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
  shape_functions(icomp,:,indices,:) = sp_col_scalar(icomp).shape_functions;
end

if (gradient || curl || divergence)
  shape_function_gradients = zeros (sp.ncomp_param, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:space.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    shape_function_gradients(icomp,:,:,indices,:) = sp_col_scalar(icomp).shape_function_gradients;
  end
end

% From here we apply the transformation
jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, msh.nel);

if (value)
  sp.shape_functions = geopdes_prod__ (msh.geo_map_jac, shape_functions);
  jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
  sp.shape_functions = bsxfun (@rdivide, sp.shape_functions, jacdet);
  jacdet = reshape (jacdet, msh.nqn, msh.nel);
end


if (divergence && ~(gradient || curl))
  sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:space.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    sp.shape_function_divs(:,indices,:) = ...
      reshape (sp_col_scalar(icomp).shape_function_gradients(icomp,:,:,:), msh.nqn, sp_col_scalar(icomp).nsh_max, msh.nel);
  end
  jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
  sp.shape_function_divs = bsxfun (@rdivide, sp.shape_function_divs, jacdet);
  jacdet = reshape (jacdet, msh.nqn, msh.nel);
end


if ((gradient || curl) && isfield (msh, 'geo_map_der2'))
%   if (msh.ndim ~= msh.rdim)
%     error ('The computation of the gradient has not been implemented for manifolds')
%   end
  
  shape_fun_grads = piola_transform_grad__ (sp, msh, shape_functions, shape_function_gradients);
  
  if (gradient)
    sp.shape_function_gradients = shape_fun_grads;
  end

  if (divergence)
    sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
    for icomp = 1:space.ncomp
      sp.shape_function_divs = sp.shape_function_divs + ...
        reshape (shape_fun_grads(icomp,icomp,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
    end
  end

  if (curl)
    if (space.ncomp == 2)
      sp.shape_function_curls = reshape (shape_fun_grads(2,1,:,:,:) - ...
	 	       shape_fun_grads(1,2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
    elseif (space.ncomp == 3)
      shape_fun_curls = zeros (space.ncomp, msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:space.ncomp
        ind1 = mod(icomp,3) + 1;
        ind2 = mod(ind1, 3) + 1;
        shape_fun_curls(icomp,:,:,:) = reshape (shape_fun_grads(ind2,ind1,:,:,:) - ...
                     shape_fun_grads(ind1,ind2,:,:,:), 1, msh.nqn, sp.nsh_max, msh.nel);
      end
      sp.shape_function_curls = shape_fun_curls;
    end
  end
end

end
