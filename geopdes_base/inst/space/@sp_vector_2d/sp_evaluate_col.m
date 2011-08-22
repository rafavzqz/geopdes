% SP_EVALUATE_COL: compute the basis functions in one column of the mesh.
%
%     sp = sp_evaluate_col (space, msh_col, 'option1', value1, ...)
%
% INPUTS:
%     
%    sp:      object defining the space of discrete functions (see sp_vector_2d)
%    msh_col: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_2d/msh_evaluate_col)
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
%    FIELD_NAME      (SIZE)                                     DESCRIPTION
%    ncomp           (scalar)                                   number of components of the functions of the space (actually, 2)
%    ndof            (scalar)                                   total number of degrees of freedom
%    ndof_dir        (2 x 2 matrix)                             for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                                   maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)                   actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)             indices of basis functions that do not vanish in each element
%    shape_functions (2 x msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%               (2 x 2 x msh_col.nqn x nsh_max x msh_col.nel)   basis function gradients evaluated at each quadrature node in each element
%    shape_function_divs (msh_col.nqn x nsh_max x msh_col.nel)  basis function gradients evaluated at each quadrature node in each element
%    shape_function_curls (msh_col.nqn x nsh_max x msh_col.nel) basis function gradients evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
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
sp1_col = sp_evaluate_col_param (space.sp1, msh, 'value', value, 'gradient', first_der);
sp2_col = sp_evaluate_col_param (space.sp2, msh, 'value', value, 'gradient', first_der);

ndof     = sp1_col.ndof + sp2_col.ndof;
ndof_dir = [sp1_col.ndof_dir; sp2_col.ndof_dir];
nsh      = sp1_col.nsh(:)' + sp2_col.nsh(:)';

connectivity = [sp1_col.connectivity; sp2_col.connectivity+sp1_col.ndof];

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
            'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
            'ncomp', 2);

% From here it depends on the transformation
if (value)
  sp.shape_functions = zeros (2, msh.nqn, sp.nsh_max, msh.nel);
  sp.shape_functions(1,:,1:sp1_col.nsh_max,:)            = sp1_col.shape_functions;
  sp.shape_functions(2,:,sp1_col.nsh_max+1:sp.nsh_max,:) = sp2_col.shape_functions;
end

if (gradient || curl || divergence)
  shape_fun_grads = zeros (2, 2, msh.nqn, sp.nsh_max, msh.nel);

  JinvT = geopdes_invT__ (msh.geo_map_jac);
  JinvT = reshape (JinvT, [2, 2, msh.nqn, msh.nel]);
  shape_fun_grads(1,:,:,1:sp1_col.nsh_max,:) = ...
              geopdes_prod__ (JinvT, sp1_col.shape_function_gradients);
  shape_fun_grads(2,:,:,sp1_col.nsh_max+1:sp.nsh_max,:) = ...
              geopdes_prod__ (JinvT, sp2_col.shape_function_gradients);

  if (gradient)
    sp.shape_function_gradients = shape_fun_grads;
  end

  if (divergence)
    sp.shape_function_divs = reshape (shape_fun_grads(1,1,:,:,:) + ...
				      shape_fun_grads(2,2,:,:,:), ...
                                      msh.nqn, sp.nsh_max, msh.nel);
  end

  if (curl)
    sp.shape_function_curls = reshape (shape_fun_grads(2,1,:,:,:) - ...
			 	       shape_fun_grads(1,2,:,:,:), ...
                                       msh.nqn, sp.nsh_max, msh.nel);
  end
end

end
