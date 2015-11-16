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
%     npatch        (scalar)                  number of patches
%     ndim          (scalar)                  dimension of the parametric space
%     rdim          (scalar)                  dimension of the physical space
%     nel           (scalar)                  number of elements in the list
%     elem_list     (1 x nel)                 numbering of the elements in the list
%     msh_patch     (1 x npatch cell-array)   evaluated elements on each patch
%
%  For more details, see the documentation
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

function msh_col = msh_evaluate_element_list (msh, elem_list, varargin)

  elem_list = elem_list(:)';

  msh_col.npatch = msh.npatch;
  msh_col.ndim = msh.ndim;
  msh_col.rdim = msh.rdim;

  msh_col.nel  = numel (elem_list);
  msh_col.elem_list = elem_list;

  if (isempty (elem_list)), return, end
  
  Nelem = cumsum ([0, msh.nel_per_patch]);
  for iptc = 1:msh.npatch
    [~,indices,~] = intersect ((Nelem(iptc)+1):Nelem(iptc+1), elem_list);
    msh_col.msh_patch{iptc} = msh_evaluate_element_list (msh.msh_patch{iptc}, indices);    
  end

end