% SP_PRECOMPUTE: precompute all the fields, as in the space structure of the technical report.
%
%     space = sp_precompute (space, msh)
%     space = sp_precompute (space, msh, 'option')
%
% INPUT:
%     
%    space: object representing the discrete function space (see sp_nurbs_2d).
%
% OUTPUT:
%
%    space: object representing the discrete function space, plus the following fields (or some of them):
%
%    FIELD_NAME      (SIZE)                                 DESCRIPTION
%    connectivity    (nsh_max x msh_col.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%             (2 x msh_col.nqn x nsh_max x msh_col.nel)     basis function gradients evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
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

function sp = sp_precompute (sp, msh, varargin)

  if (nargin == 2)
    value = true;
    gradient = true;
  else
    value = false;
    gradient = false;
    for ii=1:length(varargin)
      if (strcmpi (varargin {ii}, 'connectivity'))
        value = true;
      elseif (strcmpi (varargin {ii}, 'value'))
        value = true;
      elseif (strcmpi (varargin {ii}, 'gradient'))
        gradient = true;
      else
        error ('sp_precompute: unknown option %s', varargin {ii});
      end
    end    
  end

  if (isempty (varargin))
    sp = sp_precompute_param (sp, msh);
  else
    sp = sp_precompute_param (sp, msh, varargin{:});
  end

  if (gradient)
    if (isempty (msh.geo_map_jac))
      msh = msh_precompute (msh, 'geo_map_jac');
    end
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    JinvT = reshape (JinvT, [2, 2, msh.nqn, msh.nel]);
    sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
  end

end
