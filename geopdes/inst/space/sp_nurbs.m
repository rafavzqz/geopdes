% SP_NURBS: Constructor of the class of a tensor-product spaces of NURBS.
%
%     sp = sp_nurbs (nurbs, msh)
%     sp = sp_nurbs (knots, degree, weights, msh)
%
% INPUTS:
%
%     nurbs:     nurbs structure from the NURBS toolbox (see nrbmak)
%     msh:       msh object that defines the quadrature rule (see msh_cartesian)
%     knots:     open knot vector (cell array of size [1, ndim])
%     degree:    nurbs polynomial degree (vector of size [1, ndim])
%     weights:   weights associated to the basis functions
%
% OUTPUT:
%
%    sp: object representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        knots           (1 x ndim cell array)             knot vector in each parametric direction
%        degree          (1 x ndim vector)                 splines degree in each parametric direction
%        weights         (size(ndof_dir) matrix)           weights associated to the control points of the geometry
%        sp_univ         (1 x ndim struct)                 space of univariate splines in each parametric direction
%        ndof            (scalar)                          total number of degrees of freedom
%        ndof_dir        (1 x ndim vector)                 degrees of freedom along each direction (only for tensor product spaces)
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh_dir         (1 x ndim vector)                 maximum number of univariate shape functions per element in each parametric direction
%        boundary        (1 x 2*ndim struct array)         struct array representing the space of traces of basis functions on each edge
%
%       METHOD_NAME
%       sp_evaluate_col: compute the basis functions (and derivatives) in one column of the mesh (that is, fixing the element in the first parametric direction).
%       sp_evaluate_col_param: compute the basis functions (and derivatives) in one column of the mesh in the reference domain.
%       sp_eval_boundary_side: evaluate the basis functions in one side of the boundary.
%       sp_precompute:  compute any of the fields related to the discrete
%                       space (except boundary), in all the quadrature points,
%                       as in the space structure from previous versions.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function sp = sp_nurbs (varargin)

  if (nargin <= 3)
    nurbs = varargin{1};
    msh   = varargin{2};
    if (nargin == 3)
      transform = varargin{3};
    else
      transform = 'grad-preserving';
    end

    knots   = nurbs.knots;
    degree  = nurbs.order - 1;
    weights = squeeze (nurbs.coefs(4, :, :, :));
  else
    knots   = varargin{1};
    degree  = varargin{2};
    weights = varargin{3};
    msh     = varargin{4};
    if (nargin == 5)
      transform = varargin{5};
    else
      transform = 'grad-preserving';
    end
  end

  sp = sp_scalar (knots, degree, weights, msh, transform);

end
