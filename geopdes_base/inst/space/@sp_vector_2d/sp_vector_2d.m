% SP_VECTOR_2D: Construct the class of a two-dimensional vectorial space, using a component-wise mapping.
%
%     sp = sp_vector_2d (sp1, sp2)
%
% INPUTS:
%
%    sp1: space class of the first component (see sp_bspline_2d and sp_nurbs_2d)
%    sp2: space class of the second component (see sp_bspline_2d and sp_nurbs_2d)
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

function sp = sp_vector_2d (sp1, sp2, msh)

  sp.sp1 = sp1;
  sp.sp2 = sp2;

  sp.ncomp = 2;
  sp.nsh_max      = sp1.nsh_max + sp2.nsh_max;
  sp.ndof         = sp1.ndof + sp2.ndof;
  sp.comp_dofs{1} = 1:sp1.ndof;
  sp.comp_dofs{2} = sp1.ndof+(1:sp2.ndof);
  sp.ndof_dir     = [sp1.ndof_dir; sp2.ndof_dir];
  sp.connectivity = [sp1.connectivity; sp2.connectivity+sp1.ndof];

% For the boundary we still store everything
  if (isfield (msh, 'boundary'))
    for iside = 1:numel(msh.boundary)
      sp_bnd1 = sp1.boundary(iside);
      sp_bnd2 = sp2.boundary(iside);

      boundary.ncomp = 2;
      boundary.nsh_max      = sp_bnd1.nsh_max + sp_bnd2.nsh_max;
      boundary.nsh          = sp_bnd1.nsh + sp_bnd2.nsh;
      boundary.ndof         = sp_bnd1.ndof + sp_bnd2.ndof;
      boundary.dofs         = [sp_bnd1.dofs, sp_bnd2.dofs+sp1.ndof];

      boundary.comp_dofs{1} = sp_bnd1.dofs;
      boundary.comp_dofs{2} = sp1.ndof + sp_bnd2.dofs;
      boundary.connectivity = [sp_bnd1.connectivity; sp_bnd2.connectivity+sp_bnd1.ndof];

      boundary.shape_functions = zeros (2, msh.boundary(iside).nqn, ...
                             boundary.nsh_max, msh.boundary(iside).nel);
      boundary.shape_functions(1,:,1:sp_bnd1.nsh_max,:) = ...
                                             sp_bnd1.shape_functions;
      boundary.shape_functions(2,:, sp_bnd1.nsh_max+(1:sp_bnd2.nsh_max),:) = ...
                                             sp_bnd2.shape_functions;

      sp.boundary(iside) = boundary;
    end
  else
    sp.boundary = [];
  end

% And some function handles to rebuild the functions in different meshes
  if (isa (sp1, 'sp_bspline_2d'))
    sp.comp1 = @(MSH) sp_bspline_2d (sp1.knots, sp1.degree, MSH);
  elseif (isa (sp1, 'sp_nurbs_2d'))
    sp.comp1 = @(MSH) sp_nurbs_2d (sp1.knots, sp1.degree, sp1.weights, MSH);
  end

  if (isa (sp2, 'sp_bspline_2d'))
    sp.comp2 = @(MSH) sp_bspline_2d (sp2.knots, sp2.degree, MSH);
  elseif (isa (sp2, 'sp_nurbs_2d'))
    sp.comp2 = @(MSH) sp_nurbs_2d (sp2.knots, sp2.degree, sp2.weights, MSH);
  end

  sp = class (sp, 'sp_vector_2d');

end
