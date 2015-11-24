% SP_SCALAR: Constructor of the class of scalar tensor-product spaces (B-Splines or NURBS).
%
%     sp = sp_scalar (knots, degree, weights, msh, [transform])
%
% INPUTS:
%     
%     knots:     open knot vector (cell array of size [1, ndim])
%     degree:    spline polynomial degree (vector of size [1, ndim])
%     weights:   weights associated to the basis functions. For B-splines it should be empty
%     msh:       msh object that defines the quadrature rule (see msh_cartesian)
%     transform: string with the transformation to the physical domain, one of 
%                 'grad-preserving' (default) and 'integral-preserving', for N-forms.
%
% OUTPUT:
%
%    sp: object representing the discrete function space, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                      DESCRIPTION
%        space_type      (string)                    one of 'spline' and 'NURBS'
%        knots           (1 x ndim cell array)       knot vector in each parametric direction
%        degree          (1 x ndim vector)           splines degree in each parametric direction
%        weights         (size(ndof_dir) matrix)     weights associated to the control points of the geometry
%        sp_univ         (1 x ndim struct)           space of univariate splines in each parametric direction
%        ndof            (scalar)                    total number of degrees of freedom
%        ndof_dir        (1 x ndim vector)           number of degrees of freedom along each direction
%        nsh_max         (scalar)                    maximum number of shape functions per element
%        nsh_dir         (1 x ndim vector)           maximum number of univariate shape functions per element in each parametric direction
%        ncomp           (scalar)                    number of components of the functions of the space (actually, 1)
%        boundary        (1 x 2*ndim struct array)   struct array representing the space of traces of basis functions on each boundary side
%        transform       (string)                    one between 'grad-preserving' and 'integral-preserving'
%        dofs            (1 x ndof vector)           only for boundary spaces, degrees of freedom that do not vanish on the boundary
%        adjacent_dofs   (1 x ndof vector)           only for boundary spaces, degrees of freedom that vanish on the boundary, but its
%                                                     first derivative does not
%        constructor     function handle             function handle to construct the same discrete space in a different msh
%
%       METHODS
%       Methods that give a structure with all the functions computed in a certain subset of the mesh
%         sp_evaluate_element_list: compute basis functions (and derivatives) in a given list of elements
%         sp_evaluate_col:          compute basis functions (and derivatives) in one column of the mesh, i.e., fixing the element in the first
%                                    parametric direction
%         sp_precompute:            compute basis functions and derivatives in the whole mesh (memory consuming)       
%         sp_eval_boundary_side:    evaluate the basis functions in one side of the boundary.
%
%       Methods for post-processing, which require a computed vector of degrees of freedom
%         sp_h1_error: compute the error in H1 norm
%         sp_l2_error: compute the error in L2 norm
%         sp_eval:     evaluate the computed solution in a Cartesian grid of points
%         sp_to_vtk:   export the computed solution to a vtk file, using a Cartesian grid of points
%
%       Methods for basic connectivity operations
%         sp_get_basis_functions: compute the functions that do not vanish in a given list of elements
%         sp_get_cells:           compute the cells on which a list of functions do not vanish
%         sp_get_neighbors:       compute the neighbors, functions that share at least one element with a given one
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

function sp = sp_scalar (knots, degree, weights, msh, transform)

  if (nargin == 4)
    transform = 'grad-preserving';
  end
  
  if (isempty (weights))
    sp.space_type = 'spline';
  else
    sp.space_type = 'NURBS';
  end

% For the 1D case
  if (~iscell (knots))
    knots = {knots};
  end

  sp.knots = knots;
  sp.degree = degree;
  sp.weights = weights;

  nodes = msh.qn;
  for idim = 1:msh.ndim
    sp.sp_univ(idim) = sp_bspline_1d_param (knots{idim}, degree(idim), nodes{idim}, 'gradient', true, 'hessian', true);
  end

  sp.nsh_dir  = [sp.sp_univ.nsh_max];
  sp.nsh_max  = prod (sp.nsh_dir);
  sp.ndof_dir = [sp.sp_univ.ndof];
  sp.ndof     = prod (sp.ndof_dir);
  sp.ncomp    = 1;

% Fix a bug for 1D geometries
  if (~isempty (sp.weights))
    sp.weights = reshape (sp.weights, [sp.ndof_dir, 1]);
  end
  
% Boundary construction  
  if (msh.ndim > 1 && strcmpi (transform, 'grad-preserving'))
    for iside = 1:2*msh.ndim
%%    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
%%    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
      ind2 = ceil (iside/2);
      ind = setdiff (1:msh.ndim, ind2);

      if (~isempty (msh.boundary))
        if (strcmpi (sp.space_type, 'spline'))
          weights = [];
        elseif (strcmpi (sp.space_type, 'nurbs'))
          indices = arrayfun (@(x) 1:x, sp.ndof_dir, 'UniformOutput', false);
          if (rem (iside, 2) == 0)
            indices{ind2} = sp.ndof_dir(ind2);
          else
            indices{ind2} = 1;
          end
          weights = squeeze (sp.weights(indices{:}));
        end
        sp.boundary(iside) = sp_scalar (sp.knots(ind), sp.degree(ind), weights, msh.boundary(iside));
      end
      
      bnd_ndof_dir = sp.ndof_dir(ind);
      bnd_ndof = prod (bnd_ndof_dir);
      [ind_univ{ind}] = ind2sub (bnd_ndof_dir, 1:bnd_ndof);
      if (rem (iside, 2) == 0)
        ind_univ{ind2} = sp.ndof_dir(ind2) * ones (1, bnd_ndof);
      else
        ind_univ{ind2} = ones (1, bnd_ndof);
      end
      sp.boundary(iside).dofs = sub2ind (sp.ndof_dir, ind_univ{:});
      sp.boundary(iside).ndof = numel (sp.boundary(iside).dofs);

      if (sp.ndof_dir(ind2) > 1)
        if (rem (iside, 2) == 0)
          ind_univ{ind2} = (sp.ndof_dir(ind2) - 1) * ones (1, bnd_ndof);
        else
          ind_univ{ind2} = 2 * ones (1, bnd_ndof);
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
    sp.boundary(1).ndof = 1;
    sp.boundary(2).ndof = 1;
    
  else
    sp.boundary = [];
  end

  sp.dofs = [];
  sp.adjacent_dofs = [];
  
  sp.transform = transform;

  sp.constructor = @(MSH) sp_scalar (sp.knots, sp.degree, sp.weights, MSH, sp.transform);
  sp = class (sp, 'sp_scalar');

end
