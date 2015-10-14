% SP_PRECOMPUTE_PARAM: precompute all the fields, as in the space structure of the technical report.
%
%     space = sp_precompute_param (space, msh);
%     space = sp_precompute_param (space, msh, 'option', value);
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_bspline or sp_nurbs).
%    msh: mesh object containing the quadrature information (see msh_cartesian)
%    'option', value: additional optional parameters, available options are:
%        nsh, connectivity, value (shape_functions),
%        gradient (shape_function_gradients), divergence (shape_function_divs),
%        curl (shape_function_curls).
%     The value must be true or false. All the values are false by default.
%
% OUTPUT:
%
%    space: object containing the information of the input object, plus the 
%            fields of the old structure, that are listed below. If no option
%            is given all the fields are computed. If an option is given,
%            only the selected fields will be computed.
%
%    FIELD_NAME      (SIZE)                               DESCRIPTION
%    nsh             (1 x msh.nel vector)                  actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)            indices of basis functions that do not vanish in each element
%    shape_functions (ncomp x msh.nqn x nsh_max x msh.nel) basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%             (ncomp x rdim x msh.nqn x nsh_max x msh.nel) basis function gradients evaluated at each quadrature node in each element
%    shape_function_divs (msh.nqn x nsh_max x msh.nel)     basis function divergence evaluated at each quadrature node in each element
%    shape_function_curls
%         2D:  (msh.nqn x nsh_max x msh.nel)               basis function curl evaluated at each quadrature node in each element
%         3D:  (3 x msh.nqn x nsh_max x msh.nel)        
%
% Copyright (C) 2009, 2010 Carlo de Falco
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

function sp_out = sp_precompute_param (sp, msh, varargin)

if (isempty (varargin))
  value = true;
  gradient = false;
  divergence = false;
  curl = false;
else
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
      error ('sp_precompute: unknown option %s', varargin {ii});
    end
  end
end

if (~isstruct (msh))
  msh = msh_precompute (msh);
end
  
first_der = gradient || divergence || curl;
for icomp = 1:sp.ncomp
  sp_scalar(icomp) = sp_precompute_param (sp.scalar_spaces{icomp}, msh, 'value', value, 'gradient', first_der);
end

ndof_scalar = [sp_scalar.ndof];

ndof = sum (ndof_scalar);
nsh  = zeros (1, msh.nel);
connectivity = [];
aux = 0;
for icomp = 1:sp.ncomp
  ndof_dir(icomp,:) = sp_scalar(icomp).ndof_dir;
  nsh = nsh + sp_scalar(icomp).nsh(:)';
  
  connectivity = [connectivity; sp_scalar(icomp).connectivity+aux];
  aux = aux + ndof_scalar(icomp);
end

sp_out = struct('nsh_max', sp.nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
                'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
                'ncomp', sp.ncomp);

if (value)
  sp_out.shape_functions = zeros (sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
  for icomp = 1:sp.ncomp
    indices = sp.cumsum_nsh(icomp)+(1:sp_scalar(icomp).nsh_max);
    sp_out.shape_functions(icomp,:,indices,:) = sp_scalar(icomp).shape_functions;
  end
end

if (gradient || divergence || curl)
  shape_fun_grads = zeros (sp.ncomp, msh.rdim, msh.nqn, sp.nsh_max, msh.nel);
  aux = 0;
  for icomp = 1:sp.ncomp
    shape_fun_grads(icomp,:,:,aux+(1:sp_scalar(icomp).nsh_max),:) = ...
        scalar_sp{icomp}.shape_function_gradients;
    aux = aux + scalar_sp{icomp}.nsh_max;
  end    

  if (gradient)
    sp_out.shape_function_gradients = shape_fun_grads;
  end
    
  if (divergence)
    sp_out.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
    for icomp = 1:sp.ncomp
      sp_out.shape_function_divs = sp.shape_function_divs + ...
          reshape (shape_fun_grads(icomp,icomp,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
    end
  end
    
  if (curl)
    if (sp.ncomp == 2)
      sp_out.shape_function_curls = reshape (shape_fun_grads(2,1,:,:,:) - ...
 		 	       shape_fun_grads(1,2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
    elseif (sp.ncomp == 3)
      shape_fun_curls = zeros (sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:sp.ncomp
        ind1 = mod(icomp,3) + 1;
        ind2 = mod(ind1, 3) + 1;
        shape_fun_curls(icomp,:,:,:) = reshape (shape_fun_grads(ind2,ind1,:,:,:) - ...
                   shape_fun_grads(ind1,ind2,:,:,:), 1, msh.nqn, sp.nsh_max, msh.nel);
      end
      sp_out.shape_function_curls = shape_fun_curls;
    end
  end
end

if (isfield (struct (sp), 'dofs'))
  sp_out.dofs = sp.dofs;
end

end
