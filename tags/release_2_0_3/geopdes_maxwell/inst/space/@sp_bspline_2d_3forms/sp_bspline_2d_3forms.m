% SP_BSPLINE_2D_3FORMS: Constructor of the class of a tensor-product space of B-Splines in 2D, using a transformation for 3-forms.
%
%     sp = sp_bspline_2d_3forms (knots, degree, msh)
%
% INPUTS:
%     
%     knots:     open knot vector    
%     degree:    b-spline polynomial degree
%     msh:       msh object that defines the quadrature rule (see msh_2d)
%
% OUTPUT:
%
%    sp: object representing the discrete function space, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                      DESCRIPTION
%        spu             (struct)                    space of univariate splines in the first parametric direction
%        spv             (struct)                    space of univariate splines in the second parametric direction
%        ndof            (scalar)                    total number of degrees of freedom
%        ndof_dir        (1 x 2 vector)              degrees of freedom along each direction
%        nsh_max         (scalar)                    maximum number of shape functions per element
%        nsh_dir         (1 x 2 vector)              maximum number of univariate shape functions per element in each parametric direction
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

function sp = sp_bspline_2d_3forms (knots, degree, msh)

  sp.knots = knots;
  sp.degree = degree;

  nodes = msh.qn;
  sp.spu = sp_bspline_1d_param (knots{1}, degree(1), nodes{1}, 'gradient', false, 'hessian', false);
  sp.spv = sp_bspline_1d_param (knots{2}, degree(2), nodes{2}, 'gradient', false, 'hessian', false);
  
  sp.nsh_max  = sp.spu.nsh_max * sp.spv.nsh_max;
  sp.nsh_dir  = [sp.spu.nsh_max, sp.spv.nsh_max];
  sp.ndof     = sp.spu.ndof * sp.spv.ndof;
  sp.ndof_dir = [sp.spu.ndof, sp.spv.ndof];
  sp.ncomp    = 1;

  sp.constructor = @(MSH) sp_bspline_2d_3forms (sp.knots, sp.degree, MSH);
  sp = class (sp, 'sp_bspline_2d_3forms');

end
