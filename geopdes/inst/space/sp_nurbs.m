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
%    sp: object of the class sp_scalar, representing the discrete function space. See sp_scalar for details.
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
