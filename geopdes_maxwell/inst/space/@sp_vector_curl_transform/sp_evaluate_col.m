% SP_EVALUATE_COL: compute the basis functions in one column of the mesh.
%
%     sp = sp_evaluate_col (space, msh_col, 'option1', value1, ...)
%
% INPUTS:
%     
%    sp:      object defining the space of discrete functions (see sp_vector_curl_transform)
%    msh_col: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_cartesian/msh_evaluate_col)
%    option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            curl       |      false      |  compute shape_function_curls
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                       DESCRIPTION
%    ncomp           (scalar)                                     number of components of the functions of the space (2 or 3)
%    ndof            (scalar)                                     total number of degrees of freedom
%    ndof_dir        (ndim x ndim matrix)                         for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                                     maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)                     actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)               indices of basis functions that do not vanish in each element
%    shape_functions (rdim x msh_col.nqn x nsh_max x msh_col.nel) basis functions evaluated at each quadrature node in each element
%  In 2D (ndim = 2)
%    shape_function_curls (msh_col.nqn x nsh_max x msh_col.nel)   basis function curls evaluated at each quadrature node in each element
%  In 3D (ndim = 3)
%    shape_function_curls (rdim x msh_col.nqn x nsh_max x msh_col.nel) basis function curls evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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
curl = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'curl'))
      curl = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient') || strcmpi (varargin {ii}, 'divergence'))
%      warning ('Gradient and diveregence not implemented yet')
    else
      error ('sp_evaluate_col: unknown option %s', varargin {ii});
    end
  end
end

first_der = curl;
for icomp = 1:space.ncomp_param
  sp_col_scalar(icomp) = sp_evaluate_col_param (space.scalar_spaces{icomp}, msh, 'value', value, 'gradient', first_der);
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

% And now we apply the transformation
[JinvT, jacdet] = geopdes_invT__ (msh.geo_map_jac);

if (value)
  JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
  shape_functions = zeros (msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:sp.ncomp_param
    indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
    shape_functions(icomp,:,indices,:) = sp_col_scalar(icomp).shape_functions;
  end
  sp.shape_functions = geopdes_prod__ (JinvT, shape_functions);
end

if (curl)
  if (sp.ncomp_param == 2)
    shape_function_curls = zeros (msh.nqn, sp.nsh_max, msh.nel);
    for icomp = 1:sp.ncomp_param
      ind = setdiff (1:msh.ndim, icomp);
      indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);

      shape_function_curls(:,indices,:) = (-1)^icomp * ...
        reshape (sp_col_scalar(icomp).shape_function_gradients(ind,:,:,:), msh.nqn, sp_col_scalar(icomp).nsh_max, msh.nel);
    end
    jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
    sp.shape_function_curls = bsxfun (@rdivide, shape_function_curls, jacdet);
  
  elseif (sp.ncomp_param == 3)
    shape_function_curls = zeros (sp.ncomp_param, msh.nqn, sp.nsh_max, msh.nel);
    for icomp = 1:sp.ncomp_param
      ind(1) = mod (icomp, 3) + 1; ind(2) = mod (ind(1), 3) + 1;
      ind2 = fliplr (ind);
      indices = space.cumsum_nsh(icomp)+(1:sp_col_scalar(icomp).nsh_max);
      for ii = 1:numel(ind)
        shape_function_curls(ind(ii),:,indices,:) = (-1)^(ii-1) * ...
          reshape (sp_col_scalar(icomp).shape_function_gradients(ind2(ii),:,:,:), msh.nqn, sp_col_scalar(icomp).nsh_max, msh.nel);
      end
    end
    
    jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
    shape_function_curls = geopdes_prod__ (msh.geo_map_jac, shape_function_curls);
    sp.shape_function_curls = bsxfun (@rdivide, shape_function_curls, jacdet);
  end

end

end
