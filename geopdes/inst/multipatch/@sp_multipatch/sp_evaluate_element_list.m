% SP_EVALUATE_ELEMENT_LIST: compute the basis functions in a given list of elements.
%
%     sp = sp_evaluate_element_list (space, msh_elems, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:     object defining the space of discrete functions (see sp_multipatch)
%    msh_elems: structure containing the information of quadrature or
%               visualization points, for a given list of elements and patches
%               (see msh_multipatch/msh_evaluate_element_list)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians (only scalars)
%            laplacian  |      false      |  compute shape_function_laplacians (only scalars)
%            div        |      false      |  compute shape_function_divs (only vectors)
%            curl       |      false      |  compute shape_function_curls (only vectors)
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%
%    FIELD_NAME      (SIZE)                    DESCRIPTION
%    npatch          (scalar)                  number of patches
%    ncomp           (scalar)                  number of components of the functions of the space
%    ndof            (scalar)                  total number of degrees of freedom
%    sp_patch        (1 x npatch cell-array)   the computed space structures for each patch
%
% The information for each patch is stored in a cell-array, to allow
%  different degrees (or quadrature rules) in each direction
%
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

function sp = sp_evaluate_element_list (space, msh, varargin)

  sp.npatch = space.npatch;
  sp.ncomp = space.ncomp;
  sp.ndof = space.ndof;

  if (isempty (msh.elem_list)), return, end

  for iptc = 1:space.npatch
    sp.sp_patch{iptc} = sp_evaluate_element_list (space.sp_patch{iptc}, msh.msh_patch{iptc}, varargin{:});
    
    if (~isempty (dofs_ornt))
      error ('You have to multiply by the orientation')
    end
    sp.sp_patch{iptc}.connectivity = space.gnum{iptc}(sp.sp_patch{iptc}.connectivity);
  end

end