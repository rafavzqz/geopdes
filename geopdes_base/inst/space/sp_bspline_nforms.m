% SP_BSPLINE_NFORMS: Constructor of the class of a tensor-product space of B-Splines, using a transformation for N-forms.
%
%     sp = sp_bspline_nforms (knots, degree, msh)
%
% INPUTS:
%     
%     knots:     open knot vector (cell array of size [1, ndim])
%     degree:    b-spline polynomial degree (vector of size [1, ndim])
%     msh:       msh object that defines the quadrature rule (see msh_cartesian)
%
% OUTPUT:
%
%    sp: object representing the discrete function space, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                      DESCRIPTION
%        knots           (1 x ndim cell array)       knot vector in each parametric direction
%        degree          (1 x ndim vector)           splines degree in each parametric direction
%        sp_univ         (1 x ndim struct)           space of univariate splines in each parametric direction
%        ndof            (scalar)                    total number of degrees of freedom
%        ndof_dir        (1 x ndim vector)           number of degrees of freedom along each direction
%        nsh_max         (scalar)                    maximum number of shape functions per element
%        nsh_dir         (1 x ndim vector)           maximum number of univariate shape functions per element in each parametric direction
%        ncomp           (scalar)                    number of components of the functions of the space (actually, 1)
%
%       METHOD_NAME
%       sp_evaluate_col: compute the basis functions in one column of the mesh (that is, fixing the element in the first parametric direction).
%       sp_evaluate_col_param: compute the basis functions in one column of the mesh in the reference domain.
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

function sp = sp_bspline_nforms (knots, degree, msh)

  sp = sp_scalar (knots, degree, [], msh, 'integral-preserving');

end
