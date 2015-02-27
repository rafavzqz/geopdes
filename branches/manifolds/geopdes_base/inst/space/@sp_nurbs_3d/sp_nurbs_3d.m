% SP_NURBS_3D: Constructor of the class of a tensor-product spaces of NURBS in 3D.
%
%     sp = sp_nurbs_3d (nurbs, msh)
%     sp = sp_nurbs_3d (knots, degree, weights, msh)
%
% INPUTS:
%
%     nurbs:     nurbs structure representing a volume
%     msh:       msh object that defines the quadrature rule (see msh_3d)
%     knots:     open knot vector
%     degree:    nurbs polynomial degree (order minus one)
%     weights:   weights associated to the basis functions
%
% OUTPUT:
%
%    sp: object representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        knots           (1 x 3 cell array)                knot vector in each parametric direction
%        degree          (1 x 3 vector)                    splines degree in each parametric direction
%        weights         (size(ndof_dir) matrix)           weights associated to the control points of the geometry
%        spu             (struct)                          space of univariate splines in the first parametric direction
%        spv             (struct)                          space of univariate splines in the second parametric direction
%        spw             (struct)                          space of univariate splines in the third parametric direction
%        ndof            (scalar)                          total number of degrees of freedom
%        ndof_dir        (1 x 3 vector)                    degrees of freedom along each direction (only for tensor product spaces)
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh_dir         (1 x 3 vector)                    maximum number of univariate shape functions per element in each parametric direction
%        boundary        (1 x 6 struct array)              struct array representing the space of traces of basis functions on each edge
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

function sp = sp_nurbs_3d (varargin)

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
    error ('sp_nurbs_3d: wrong input arguments. See the help for usage');
  end

  nodes = msh.qn;
  sp.spu = sp_bspline_1d_param (sp.knots{1}, sp.degree(1), nodes{1}, 'gradient', true, 'hessian', true);
  sp.spv = sp_bspline_1d_param (sp.knots{2}, sp.degree(2), nodes{2}, 'gradient', true, 'hessian', true);
  sp.spw = sp_bspline_1d_param (sp.knots{3}, sp.degree(3), nodes{3}, 'gradient', true, 'hessian', true);
  
  sp.nsh_max  = sp.spu.nsh_max * sp.spv.nsh_max * sp.spw.nsh_max;
  sp.nsh_dir  = [sp.spu.nsh_max, sp.spv.nsh_max, sp.spw.nsh_max];
  sp.ndof     = sp.spu.ndof * sp.spv.ndof * sp.spw.ndof;
  sp.ndof_dir = [sp.spu.ndof, sp.spv.ndof, sp.spw.ndof];
  sp.ncomp    = 1;

  if (~isempty (msh.boundary))
    msh_bnd = msh.boundary;
    for iside = 1:numel(msh_bnd)
      ind = setdiff (1:3, ceil(iside/2)); %ind=[2 3; 2 3; 1 3; 1 3; 1 2; 1 2]
      ind2 = floor ((iside+1)/2); % ind2 = [1 1 2 2 3 3];

      boundary.ndof_dir = sp.ndof_dir(ind);
      boundary.ndof = prod (boundary.ndof_dir);
      boundary.nsh_dir = prod (sp.nsh_dir(ind));
      boundary.nsh_max = prod (boundary.nsh_dir);
      boundary.ncomp = 1;

      idx = ones (3, boundary.ndof);
      [idx(ind(1),:), idx(ind(2),:)] = ind2sub ([boundary.ndof_dir], 1:boundary.ndof);
      if (rem (iside, 2) == 0)
        idx(ind2,:) = sp.ndof_dir(ind2);
      end

      boundary.dofs = sub2ind (sp.ndof_dir, idx(1,:), idx(2,:), idx(3,:));
      sp.boundary(iside) = boundary;
    end 
  else
    sp.boundary = [];
  end

  sp.nsh = [];
  sp.connectivity = [];
  sp.shape_functions = [];
  sp.shape_function_gradients = [];

  sp.constructor = @(MSH) sp_nurbs_3d (sp.knots, sp.degree, sp.weights, MSH);
  sp = class (sp, 'sp_nurbs_3d');

end
