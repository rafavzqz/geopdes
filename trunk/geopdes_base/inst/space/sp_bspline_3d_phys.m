% SP_BSPLINE_3D_PHYS: Construct a tensor-product space of B-Splines on the physical domain in 3D.
%
%     sp = sp_bspline_3d_phys (knots, degree, msh, 'option1', value1, ...)
%
% INPUTS:
%
%     knots:  open knot vector
%     degree: b-spline polynomial degree
%     msh:    structure containing the domain partition and the quadrature rule (see msh_push_forward_3d)
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
%        spfun           (function handle)                 function to evaluate an element of the discrete function space, given the Fourier 
%                                                          coefficients and a set of points in the parametric space
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

function sp = sp_bspline_3d_phys (knots, degree, msh, varargin)

  sp = sp_bspline_3d_param (knots, degree, msh, varargin{:});

  if (isfield (sp, 'shape_function_gradients'))
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
  end

  sp.spfun  = @(MSH) sp_bspline_3d_phys (knots, degree, MSH, varargin{:});
 
end

