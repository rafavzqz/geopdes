% SP_BSPLINE_3D: Construct the class of a tensor-product space of B-Splines in 3D.
%
%     sp = sp_bspline_3d (knots, degree, msh)
%
% INPUTS:
%     
%     knots:     open knot vector    
%     degree:    b-spline polynomial degree
%     msh:       msh class containing (in the field msh.qn) the points 
%                along each parametric direction in the parametric 
%                domain at which to evaluate, i.e. quadrature points 
%                or points for visualization (see msh_3d)
%
% OUTPUT:
%
%    sp: class representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        spu             (struct)                          space of univariate splines in the first parametric direction
%        spv             (struct)                          space of univariate splines in the second parametric direction
%        spw             (struct)                          space of univariate splines in the third parametric direction
%        ndof            (scalar)                          total number of degrees of freedom
%        ndof_dir        (1 x 3 vector)                    degrees of freedom along each direction (only for tensor product spaces)
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)              actual number of shape functions per each element
%        ncomp           (scalar)                    number of components of the functions of the space (actually, 1)
%        boundary        (1 x 6 struct array)              struct array representing the space of traces of basis functions on each edge
%
%
%       METHOD_NAME
%       sp_evaluate_col: compute the basis functions in one column of the mesh (that is, fixing the element in the first parametric direction).
%       sp_evaluate_row: compute the basis functions in one row of the mesh (that is, fixing the element in the last parametric direction).
%
%       sp_eval:        evaluate a function, given by its dofs, at a given set of points.
%       sp_to_vtk:      export a function, given by its dofs, in the vtk format.
%       sp_h1_error:    evaluate the error in H^1 norm.
%       sp_l2_error:    evaluate the error in L^2 norm.
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

function sp = sp_bspline_3d (knots, degree, msh)

  sp.knots = knots;
  sp.degree = degree;

  nodes = msh.qn;
  sp.spu = sp_bspline_1d_param (knots{1}, degree(1), nodes{1}, 'gradient', true, 'hessian', true);
  sp.spv = sp_bspline_1d_param (knots{2}, degree(2), nodes{2}, 'gradient', true, 'hessian', true);
  sp.spw = sp_bspline_1d_param (knots{3}, degree(3), nodes{3}, 'gradient', true, 'hessian', true);
  
  sp.nsh_max  = sp.spu.nsh_max * sp.spv.nsh_max * sp.spw.nsh_max;
  sp.ndof     = sp.spu.ndof * sp.spv.ndof * sp.spw.ndof;
  sp.ndof_dir = [sp.spu.ndof, sp.spv.ndof, sp.spw.ndof];
  sp.ncomp    = 1;

  ucp = sp.ndof_dir(1);
  vcp = sp.ndof_dir(2); 
  wcp = sp.ndof_dir(3);
  if (~isempty (msh.boundary))
    msh_bnd = msh.boundary;
    for iside = 1:numel(msh_bnd)
      switch (iside)
        case {1}
          ind = [2, 3];
          [vidx, widx] = ind2sub ([vcp, wcp], 1:vcp*wcp);
          uidx = ones (size (vidx));
        case {2}
          ind = [2, 3];
          [vidx, widx] = ind2sub ([vcp, wcp], 1:vcp*wcp);
          uidx = ucp * ones (size (vidx));
        case {3}
          ind = [1, 3];
          [uidx, widx] = ind2sub ([ucp, wcp], 1:ucp*wcp);
          vidx = ones (size (uidx));
        case {4}
          ind = [1, 3];
          [uidx, widx] = ind2sub ([ucp, wcp], 1:ucp*wcp);
          vidx = vcp * ones (size (uidx));
        case {5}
          ind = [1, 2];
          [uidx, vidx] = ind2sub ([ucp, vcp], 1:ucp*vcp);
          widx = ones (size (vidx));
        case {6}
          ind = [1, 2];
          [uidx, vidx] = ind2sub ([ucp, vcp], 1:ucp*vcp);
          widx = wcp * ones (size (vidx));
      end
      bnd_iside = ...
        sp_bspline_2d_param (knots(ind), degree(ind), msh_bnd(iside));

      boundary      = rmfield (bnd_iside, 'shape_function_gradients');
      boundary.dofs = sub2ind ([ucp, vcp, wcp], uidx, vidx, widx);

      sp.boundary(iside) = boundary;
    end 
  else
    sp.boundary = [];
  end

  sp.constructor = @(MSH) sp_bspline_3d (sp.knots, sp.degree, MSH);
  sp = class (sp, 'sp_bspline_3d');

end
