% SP_NURBS_3D_PARAM: Construct a tensor-product space of NURBS on the parametric domain in 3D.
%                      This function is not usually meant to be invoked directly by the user but rather
%                      through sp_nurbs_3d_phys.
%
%     sp = sp_nurbs_3d_param (nrb, msh, 'option1', value1, ...)
%
% INPUTS:
%     
%     nrb:  NURBS volume structure, as in the NURBS toolbox
%     msh:  msh structure (see msh_3h_handle)
%    'option', value: additional optional parameters, currently available options are:
%
%              Name     |   Default value |  Meaning
%           ------------+-----------------+-----------
%            gradient   |      true       |  compute shape_function_gradients
%
% OUTPUT:
%
%    sp: structure representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        ndof            (scalar)                          total number of degrees of freedom
%        ndof_dir        (1 x 3 vector)                    degrees of freedom along each direction (only for tensor product spaces)
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)              actual number of shape functions per each element
%        connectivity    (nsh_max x msh.nel vector)        indices of basis functions that do not vanish in each element
%        shape_functions (msh.nqn x nsh_max x msh.nel)     basis functions evaluated at each quadrature node in each element
%        shape_function_gradients
%                        (3 x msh.nqn x nsh_max x msh.nel) basis function gradients evaluated at each quadrature node in each element
%        boundary        (1 x 6 struct array)              struct array representing the space of traces of basis functions on each face
%
%   For more details, see the documentation
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function sp = sp_nurbs_3d_param (nrb, msh, varargin)

  w  = squeeze (nrb.coefs(4,:,:));

  sp = sp_bspline_3d_param (nrb.knots, nrb.order - 1, msh, varargin{:});
  sp = bsp_2_nrb_3d__ (sp, msh, w);
  if (isfield (msh, 'boundary'))
    for iside = 1:numel(msh.boundary)
      sp.boundary(iside) = bsp_2_nrb_2d__ (sp.boundary(iside), msh.boundary(iside), w);
    end
  end

end

