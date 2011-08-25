% SP_VECTOR_3D: Constructor of the class of three-dimensional vectorial spaces with a component-wise mapping.
%
%     sp = sp_vector_3d (sp1, sp2, sp3, msh)
%
% INPUTS:
%
%    sp1: space object of the first component (see sp_bspline_3d or sp_nurbs_3d)
%    sp2: space object of the second component (see sp_bspline_3d or sp_nurbs_3d)
%    sp3: space object of the third component (see sp_bspline_3d or sp_nurbs_3d)
%    msh: mesh object that defines the domain partition and the quadrature rule (see msh_3d)
%
% OUTPUT:
%
%    sp: object representing the discrete function space of vector-valued functions, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                      DESCRIPTION
%        sp1             (space object)              space object for the first component
%        sp2             (space object)              space object for the second component
%        sp3             (space object)              space object for the third component
%        ndof            (scalar)                    total number of degrees of freedom
%        ndof_dir        (3 x 3 matrix)              for each component, number of degrees of freedom along each direction
%        comp_dofs       (1 x 3 cell array)          indices of the degrees of freedom for each component
%        nsh_max         (scalar)                    maximum number of shape functions per element
%        ncomp           (scalar)                    number of components of the functions of the space (actually, 3)
%        boundary        (1 x 6 struct array)        struct array representing the space of traces of basis functions on each edge
%
%       METHOD_NAME
%       sp_evaluate_col: compute the basis functions in one column of the mesh (that is, fixing the element in the first parametric direction).
%       sp_eval_boundary_side: evaluate the basis functions in one side of the boundary.
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

function sp = sp_vector_3d (sp1, sp2, sp3, msh)

  sp.sp1 = sp1;
  sp.sp2 = sp2;
  sp.sp3 = sp3;

  sp.ncomp = 3;
  sp.nsh_max      = sp1.nsh_max + sp2.nsh_max + sp3.nsh_max;
  sp.ndof         = sp1.ndof + sp2.ndof + sp3.ndof;
  sp.comp_dofs{1} = 1:sp1.ndof;
  sp.comp_dofs{2} = sp1.ndof + (1:sp2.ndof);
  sp.comp_dofs{3} = sp1.ndof + sp2.ndof + (1:sp3.ndof);
  sp.ndof_dir     = [sp1.ndof_dir; sp2.ndof_dir; sp3.ndof_dir];

  if (~isempty (msh.boundary))
    for iside = 1:numel(msh.boundary)
      sp_bnd1 = sp1.boundary(iside);
      sp_bnd2 = sp2.boundary(iside);
      sp_bnd3 = sp3.boundary(iside);

      boundary.ncomp = 3;
      boundary.nsh_max      = sp_bnd1.nsh_max + sp_bnd2.nsh_max + sp_bnd3.nsh_max;
      boundary.ndof         = sp_bnd1.ndof + sp_bnd2.ndof + sp_bnd3.ndof;
      boundary.ndof_dir     = [sp_bnd1.ndof_dir; sp_bnd2.ndof_dir; sp_bnd3.ndof_dir];
      boundary.dofs         = [sp_bnd1.dofs, sp_bnd2.dofs + sp1.ndof, ...
			       sp_bnd3.dofs + sp1.ndof + sp2.ndof];

      boundary.comp_dofs{1} = sp_bnd1.dofs;
      boundary.comp_dofs{2} = sp1.ndof + sp_bnd2.dofs;
      boundary.comp_dofs{3} = sp1.ndof + sp2.ndof + sp_bnd3.dofs;

      sp.boundary(iside) = boundary;
    end
  else
    sp.boundary = [];
  end

  sp.nsh = [];
  sp.connectivity = [];
  sp.shape_functions = [];
  sp.shape_function_gradients = [];
  sp.shape_function_divs = [];
  sp.shape_function_curls = [];

  sp.constructor = @(MSH) sp_vector_3d (sp1.constructor (MSH), ...
                                        sp2.constructor (MSH), ...
                                        sp3.constructor (MSH), MSH);
  sp = class (sp, 'sp_vector_3d');

end
