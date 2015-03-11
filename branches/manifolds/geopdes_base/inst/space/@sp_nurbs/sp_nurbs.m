% SP_NURBS: Constructor of the class of a tensor-product spaces of NURBS.
%
%     sp = sp_nurbs (nurbs, msh)
%     sp = sp_nurbs (knots, degree, weights, msh)
%
% INPUTS:
%
%     nurbs:     nurbs structure representing a volume
%     msh:       msh object that defines the quadrature rule (see msh_cartesian)
%     knots:     open knot vector
%     degree:    nurbs polynomial degree (order minus one)
%     weights:   weights associated to the basis functions
%
% OUTPUT:
%
%    sp: object representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        knots           (1 x ndim cell array)             knot vector in each parametric direction
%        degree          (1 x ndim vector)                 splines degree in each parametric direction
%        weights         (size(ndof_dir) matrix)           weights associated to the control points of the geometry
%        sp_univ         (1 x ndim struct)                 space of univariate splines in each parametric direction
%        ndof            (scalar)                          total number of degrees of freedom
%        ndof_dir        (1 x ndim vector)                 degrees of freedom along each direction (only for tensor product spaces)
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh_dir         (1 x ndim vector)                 maximum number of univariate shape functions per element in each parametric direction
%        boundary        (1 x 2*ndim struct array)         struct array representing the space of traces of basis functions on each edge
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

function sp = sp_nurbs (varargin)

  if (nargin == 2)
    nurbs = varargin{1};
    msh   = varargin{2};

    sp.knots = nurbs.knots;
    sp.degree = nurbs.order - 1;
    sp.weights = squeeze (nurbs.coefs(4, :, :, :));
  elseif (nargin == 4)
    sp.knots   = varargin{1};
    sp.degree  = varargin{2};
    sp.weights = varargin{3};
    msh        = varargin{4};
  else
    error ('sp_nurbs: wrong input arguments. See the help for usage');
  end

  nodes = msh.qn;
  for idim = 1:msh.ndim
    sp.sp_univ(idim) = sp_bspline_1d_param (sp.knots{idim}, sp.degree(idim), nodes{idim}, 'gradient', true, 'hessian', true);
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

      boundary(iside).ndof_dir = sp.ndof_dir(ind);
      boundary(iside).ndof = prod (boundary(iside).ndof_dir);
      boundary(iside).nsh_dir = sp.nsh_dir(ind);
      boundary(iside).nsh_max = prod (boundary(iside).nsh_dir);
      boundary(iside).ncomp = 1;
      
      [ind_univ{ind}] = ind2sub ([boundary(iside).ndof_dir], 1:boundary(iside).ndof);
      if (rem (iside, 2) == 0)
        ind_univ{ind2} = sp.ndof_dir(ind2) * ones (1, boundary(iside).ndof);
      else
        ind_univ{ind2} = ones (1, boundary(iside).ndof);
      end

      boundary(iside).dofs = sub2ind (sp.ndof_dir, ind_univ{:});

      if (sp.ndof_dir(ind2) > 1)
        if (rem (iside, 2) == 0)
          ind_univ{ind2} = (sp.ndof_dir(ind2) - 1) * ones (1, boundary(iside).ndof);
        else
          ind_univ{ind2} = 2 * ones (1, boundary(iside).ndof);
        end
        boundary(iside).adjacent_dofs = sub2ind (sp.ndof_dir, ind_univ{:});
      end
    end
        
    sp.boundary = boundary;
  else
    sp.boundary = [];
  end
  
  sp.nsh = [];
  sp.connectivity = [];
  sp.shape_functions = [];
  sp.shape_function_gradients = [];

  sp.spline_space = sp_bspline (sp.knots, sp.degree, msh);
  

  sp.constructor = @(MSH) sp_nurbs (sp.knots, sp.degree, sp.weights, MSH);
  sp = class (sp, 'sp_nurbs');
  
end
