% SP_BSPLINE: Constructor of the class of tensor-product spaces of B-Splines.
%
%     sp = sp_bspline (knots, degree, msh, [transform])
%
% INPUTS:
%     
%     knots:     open knot vector (cell array of size [1, ndim])
%     degree:    b-spline polynomial degree (vector of size [1, ndim])
%     msh:       msh object that defines the quadrature rule (see msh_cartesian)
%     transform: string with the transform to the physical domain, one of 
%                 'grad-preserving' (default) and 'integral-preserving', for N-forms.
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

function sp = sp_bspline (knots, degree, msh, transform)

if (nargin == 3)
  transform = 'grad-preserving';
end

sp = sp_scalar (knots, degree, [], msh, transform);

end