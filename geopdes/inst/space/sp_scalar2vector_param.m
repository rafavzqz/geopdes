% SP_SCALAR2VECTOR_PARAM: from a scalar space struct, compute a vector-valued
%  space struct, with the same basis functions on each component of the parametric domain.
%
%   space_vec = sp_scalar2vector_param (space, msh, 'option1', value1, ...);
%
% INPUT:
%     
%     space: struct defining the scalar discrete space (see for instance sp_scalar/sp_evaluate_col)
%     msh:   struct defining the mesh (see for instance msh_cartesian/msh_evaluate_col)
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%
% OUTPUT:
%
%     space_vec: struct defining the vector-valued discrete space, with
%                 ndof = space.ndof*msh.ndim.
% 
% Copyright (C) 2022 Rafael Vazquez
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

function space_vec = sp_scalar2vector_param (space, msh, varargin)

value = true;
gradient = false;
hessian = false;
  if (~isempty (varargin))
    if (~rem (length (varargin), 2) == 0)
      error ('sp_scalar2vector: options must be passed in the [option, value] format');
    end
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin {ii}, 'value'))
        value = varargin {ii+1};
      elseif (strcmpi (varargin {ii}, 'gradient'))
        gradient = varargin {ii+1};
      elseif (strcmpi (varargin {ii}, 'hessian'))
        hessian = varargin {ii+1};
      else
        warning ('Ignoring unknown option %s', varargin {ii});
      end
    end
  end
  
  ncomp = msh.rdim;
  space_vec.nsh_max = ncomp * space.nsh_max;
  space_vec.nsh = ncomp * space.nsh;
  space_vec.ndof = ncomp * space.ndof;
  space_vec.ndof_dir = repmat (space.ndof_dir, ncomp, 1);
  space_vec.ncomp = ncomp;

  space_vec.connectivity = [];
  for icomp = 1:ncomp
    space_vec.connectivity = [space_vec.connectivity; space.connectivity+(icomp-1)*space.ndof];
  end

  if (value)
    space_vec.shape_functions = zeros (ncomp, msh.nqn, space_vec.nsh_max, msh.nel);
    for icomp = 1:ncomp
      fun_inds = (icomp-1)*space.nsh_max+1:icomp*space.nsh_max;
      space_vec.shape_functions(icomp,:,fun_inds,:) = space.shape_functions;
    end
  end
  
  if (gradient)
    space_vec.shape_function_gradients = zeros (ncomp, msh.ndim, msh.nqn, space_vec.nsh_max, msh.nel);
    for icomp = 1:ncomp
      fun_inds = (icomp-1)*space.nsh_max+1:icomp*space.nsh_max;
      space_vec.shape_function_gradients(icomp,:,:,fun_inds,:) = space.shape_function_gradients;
    end
  end

  if (hessian)
    space_vec.shape_function_hessians = zeros (ncomp, msh.ndim, msh.ndim, msh.nqn, space_vec.nsh_max, msh.nel);
    for icomp = 1:ncomp
      fun_inds = (icomp-1)*space.nsh_max+1:icomp*space.nsh_max;
      space_vec.shape_function_hessians(icomp,:,:,:,fun_inds,:) = space.shape_function_hessians;
    end
  end
  
end
