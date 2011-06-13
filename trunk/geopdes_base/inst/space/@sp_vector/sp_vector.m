% SP_VECTOR: Construct the class of general vectorial space.
%
%     sp = sp_vector (sp1, sp2)
%     sp = sp_vector (sp1, sp2, sp3)
%
% INPUTS:
%
%    sp1: space class of the first component (see sp_bspline_2d and sp_nurbs_2d)
%    sp2: space class of the second component (see sp_bspline_2d and sp_nurbs_2d%    sp3: space class of the third component (see sp_bspline_3d and sp_nurbs_3d)
%
% OUTPUT:
%
%
%
%
%
%    sp: class representing the discrete function space, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                      DESCRIPTION
%        spu             (struct)                    space of univariate splines in the first parametric direction
%        spv             (struct)                    space of univariate splines in the second parametric direction
%        ndof            (scalar)                    total number of degrees of freedom
%        ndof_dir        (1 x 2 vector)              degrees of freedom along each direction
%        nsh_max         (scalar)                    maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)        actual number of shape functions per each element
%        connectivity    (nsh_max x msh.nel vector)  indices of basis functions that do not vanish in each element
%        ncomp           (scalar)                    number of components of the functions of the space (actually, 1)
%        boundary        (1 x 4 struct array)        struct array representing the space of traces of basis functions on each edge
%
%       METHOD_NAME
%       sp_evaluate_col: compute the basis functions in one column of the mesh (that is, fixing the element in the first parametric direction).
%
%       sp_evaluate_row: compute the basis functions in one row of the mesh (that is, fixing the element in the last parametric direction).
%
%       sp_eval:        evaluate a function, given by its dofs, at a given set of points.
%       sp_to_vtk:      export a function, given by its dofs, in the vtk format.
%       sp_h1_error:    evaluate the error in H^1 norm.
%       sp_l2_error:    evaluate the error in L^2 norm.
%       op_u_v:         compute the mass matrix.
%       op_gradu_gradv: compute the stifness matrix.
%       op_f_v:         compute the right-hand side.
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

function sp = sp_vector (varargin)

  if (nargin == 2) % Two-dimensional space
    sp.sp1 = varargin{1};
    sp.sp2 = varargin{2};

    sp.ncomp = 2;
    sp.ndof         = sp.sp1.ndof + sp.sp2.ndof;
    sp.comp_dofs{1} = 1:sp.sp1.ndof;
    sp.comp_dofs{2} = sp.sp1.ndof+(1:sp.sp2.ndof);
    sp.ndof_dir     = [sp.sp1.ndof_dir; sp.sp2.ndof_dir];

  elseif (nargin == 3) % Three-dimensional space
    sp.sp1 = varargin{1};
    sp.sp2 = varargin{2};
    sp.sp3 = varargin{3};

    sp.ncomp = 3;
    sp.ndof         = sp.sp1.ndof + sp.sp2.ndof + sp.sp3.ndof;
    sp.comp_dofs{1} = 1:sp.sp1.ndof;
    sp.comp_dofs{2} = sp.sp1.ndof + (1:sp.sp2.ndof);
    sp.comp_dofs{3} = sp.sp1.ndof + sp.sp2.ndof + (1:sp.sp3.ndof);
    sp.ndof_dir     = [sp.sp1.ndof_dir; sp.sp2.ndof_dir; sp.sp3.ndof_dir];
  end

  sp = class (sp, 'sp_vector');

end
