% SP_PRECOMPUTE: precompute all the fields, as in the space structure of the technical report.
%
%     space = sp_precompute (space, msh);
%     space = sp_precompute (space, msh, 'option', value);
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_bspline or sp_nurbs).
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

function sp = sp_precompute (sp, msh, varargin)

  if (isempty (varargin))
    nsh = true;
    connectivity = true;
    value = true;
    gradient = true;
    divergence = true;
    curl = false;
    scalar_gradient = true;
  else
    if (~rem (length (varargin), 2) == 0)
      error ('sp_precompute: options must be passed in the [option, value] format');
    end

    nsh = false;
    connectivity = false;
    value = false;
    gradient = false;
    divergence = false;
    curl = false;
    scalar_gradient = false;
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin{ii}, 'connectivity'))
        connectivity = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'nsh'))
        nsh = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'value'))
        value = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'gradient'))
        gradient = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'divergence'))
        divergence = varargin{ii+1};
      elseif (strcmpi (varargin{ii}, 'curl'))
        curl = varargin{ii+1};
      else
        error ('sp_precompute: unknown option %s', varargin {ii});
      end
    end
    if (gradient || divergence || curl)
      scalar_gradient = true;
    end
  end
  
% % Derivatives are not computed for boundary spaces (at least for now)
%   if (~isempty (sp.dofs))
%      gradient = false;
%      divergence = false;
%      curl = false;
%   end

% Precompute everything for each component (already in the physical domain)
% These structures will not be stored in memory
  scalar_sp = cell (sp.ncomp,1);
  for icomp = 1:sp.ncomp
    scalar_sp{icomp} = sp_precompute (sp.scalar_spaces{icomp}, msh, 'nsh', nsh, 'connectivity', ...
                    connectivity, 'value', value, 'gradient', scalar_gradient);
  end

  if (nsh)
    sp.nsh  = zeros (1, msh.nel);
    for icomp = 1:sp.ncomp
      sp.nsh = sp.nsh + scalar_sp{icomp}.nsh(:)';
    end
  end

  if (connectivity)
    sp.connectivity = [];
    aux = 0;
    for icomp = 1:sp.ncomp
      sp.connectivity = [sp.connectivity; scalar_sp{icomp}.connectivity+aux];
      aux = aux + scalar_sp{icomp}.ndof;
    end
  end

  if (value)
    sp.shape_functions = zeros (sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
    aux = 0;
    for icomp = 1:sp.ncomp
      sp.shape_functions(icomp,:,aux+(1:scalar_sp{icomp}.nsh_max),:) = scalar_sp{icomp}.shape_functions;
      aux = aux + scalar_sp{icomp}.nsh_max;
    end
  end

  if (gradient || divergence || curl)
    shape_fun_grads = zeros (sp.ncomp, msh.rdim, msh.nqn, sp.nsh_max, msh.nel);
    aux = 0;
    for icomp = 1:sp.ncomp
      shape_fun_grads(icomp,:,:,aux+(1:scalar_sp{icomp}.nsh_max),:) = ...
          scalar_sp{icomp}.shape_function_gradients;
      aux = aux + scalar_sp{icomp}.nsh_max;
    end    
    
    if (gradient)
      sp.shape_function_gradients = shape_fun_grads;
    end
    
    if (divergence)
      sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:sp.ncomp
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
