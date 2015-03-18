% SP_BSPLINE: Constructor of the class of tensor-product spaces of B-Splines.
%
%     sp = sp_bspline (knots, degree, msh)
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
%        boundary        (1 x 2*ndim struct array)   struct array representing the space of traces of basis functions on each edge
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

function sp = sp_bspline (knots, degree, msh)

% For the 1D case
  if (~iscell (knots))
    knots = {knots};
  end

  sp.knots = knots;
  sp.degree = degree;

  nodes = msh.qn;
  for idim = 1:msh.ndim
    sp.sp_univ(idim) = sp_bspline_1d_param (knots{idim}, degree(idim), nodes{idim}, 'gradient', true, 'hessian', true);
  end

  sp.nsh_dir  = [sp.sp_univ.nsh_max];
  sp.nsh_max  = prod (sp.nsh_dir);
  sp.ndof_dir = [sp.sp_univ.ndof];
  sp.ndof     = prod (sp.ndof_dir);
  sp.ncomp    = 1;

  if (~isempty (msh.boundary))
    for iside = 1:numel (msh.boundary)
%%    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
%%    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
      ind2 = ceil (iside/2);
      ind = setdiff (1:msh.ndim, ind2);

      sp.boundary(iside) = sp_bspline (sp.knots(ind), sp.degree(ind), msh.boundary(iside));
      
      [ind_univ{ind}] = ind2sub ([sp.boundary(iside).ndof_dir], 1:sp.boundary(iside).ndof);
      if (rem (iside, 2) == 0)
        ind_univ{ind2} = sp.ndof_dir(ind2) * ones (1, sp.boundary(iside).ndof);
      else
        ind_univ{ind2} = ones (1, sp.boundary(iside).ndof);
      end
      sp.boundary(iside).dofs = sub2ind (sp.ndof_dir, ind_univ{:});

      if (sp.ndof_dir(ind2) > 1)
        if (rem (iside, 2) == 0)
          ind_univ{ind2} = (sp.ndof_dir(ind2) - 1) * ones (1, sp.boundary(iside).ndof);
        else
          ind_univ{ind2} = 2 * ones (1, sp.boundary(iside).ndof);
        end
        sp.boundary(iside).adjacent_dofs = sub2ind (sp.ndof_dir, ind_univ{:});
      end
      
    end
        
  elseif (msh.ndim == 1)
    sp.boundary(1).dofs = 1;
    sp.boundary(2).dofs = sp.ndof;
    if (sp.ndof > 1)
      sp.boundary(1).adjacent_dofs = 2;
      sp.boundary(2).adjacent_dofs = sp.ndof - 1;
    end
    
  else
    sp.boundary = [];
  end

  sp.nsh = [];
  sp.connectivity = [];
  sp.shape_functions = [];
  sp.shape_function_gradients = [];
  sp.shape_function_hessians = [];
  sp.dofs = [];
  sp.adjacent_dofs = [];

  sp.constructor = @(MSH) sp_bspline (sp.knots, sp.degree, MSH);
  sp = class (sp, 'sp_bspline');

end
