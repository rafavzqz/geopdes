% SP_EVALUATE_COL: compute the basis functions in one column of the mesh.
%
%     sp = sp_evaluate_col (space, msh, colnum, 'option1', value1, ...)
%
% INPUTS:
%     
%     space:  class defining the space of discrete functions (see sp_nurbs_3d)
%     msh:    msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_3d/msh_evaluate_col)
%     colnum: number of the fixed element in the first parametric direction
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      true       |  compute shape_function_gradients
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%
%    FIELD_NAME      (SIZE)                            DESCRIPTION
%    ncomp           (scalar)                          number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                          total number of degrees of freedom
%    ndof_dir        (1 x 3 vector)                    degrees of freedom along each direction
%    nsh_max         (scalar)                          maximum number of shape functions per element
%    nsh             (1 x msh.nelcol vector)           actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nelcol vector)     indices of basis functions that do not vanish in each element
%    shape_functions (msh.nqn x nsh_max x msh.nelcol)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%                 (3 x msh.nqn x nsh_max x msh.nelcol) basis function gradients evaluated at each quadrature node in each element
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

function [sp, elem_list] = sp_evaluate_col (space, msh, colnum, varargin)

value = true;
gradient = true;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    else
      error ('sp_evaluate_col: unknown option %s', varargin {ii});
    end
  end
end

sp = sp_onecol_param (space, msh, colnum, varargin{:});

nel_col = msh.nelv * msh.nelw;
indu = colnum * ones(msh.nelv, msh.nelw);
indv = repmat ((1:msh.nelv)', 1, msh.nelw);
indw = repmat ((1:msh.nelw), msh.nelv, 1);

elem_list = sub2ind ([msh.nelu, msh.nelv, msh.nelw], indu, indv, indw);
elem_list = elem_list(:);

if (gradient)
  JinvT = geopdes_invT__ (msh.geo_map_jac(:,:,:,elem_list));
  JinvT = reshape (JinvT, [3, 3, msh.nqn, nel_col]);
  sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
end

end
