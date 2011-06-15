% SP_NURBS_3D: Construct the class of a tensor-product space of NURBS in 3D.
%
%     sp = sp_nurbs_3d (nurbs, msh)
%     sp = sp_nurbs_3d (knots, degree, weights, msh)
%
% INPUTS:
%
%     nurbs:     nurbs structure representing a volume
%     msh:       msh class containing (in the field msh.qn) the points 
%                along each parametric direction in the parametric 
%                domain at which to evaluate, i.e. quadrature points 
%                or points for visualization (see msh_3d)
%     knots:     open knot vector
%     degree:    nurbs polynomial degree (order minus one)
%     weights:   weights associated to the basis functions
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
%        connectivity    (nsh_max x msh.nel vector)        indices of basis functions that do not vanish in each element
%        boundary        (1 x 6 struct array)              struct array representing the space of traces of basis functions on each edge
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
  sp.ndof     = sp.spu.ndof * sp.spv.ndof * sp.spw.ndof;
  sp.ndof_dir = [sp.spu.ndof, sp.spv.ndof, sp.spw.ndof];
  sp.ncomp    = 1;

  conn_u = reshape (sp.spu.connectivity, sp.spu.nsh_max, 1, 1, msh.nelu, 1, 1);
  conn_u = repmat  (conn_u, [1, sp.spv.nsh_max, sp.spw.nsh_max, 1, msh.nelv, msh.nelw]);
  conn_u = reshape (conn_u, [], msh.nel);

  conn_v = reshape (sp.spv.connectivity, 1, sp.spv.nsh_max, 1, 1, msh.nelv, 1);
  conn_v = repmat  (conn_v, [sp.spu.nsh_max, 1, sp.spw.nsh_max, msh.nelu, 1, msh.nelw]);
  conn_v = reshape (conn_v, [], msh.nel);

  conn_w = reshape (sp.spw.connectivity, 1, 1, sp.spw.nsh_max, 1, 1, msh.nelw);
  conn_w = repmat  (conn_w, [sp.spu.nsh_max, sp.spv.nsh_max, 1, msh.nelu, msh.nelv, 1]);
  conn_w = reshape (conn_w, [], msh.nel);

  connectivity = zeros (sp.nsh_max, msh.nel);
  indices = (conn_u ~= 0) & (conn_v ~= 0);
  connectivity(indices) = ...
     sub2ind ([sp.spu.ndof, sp.spv.ndof, sp.spw.ndof], conn_u(indices), conn_v(indices), conn_w(indices));
  sp.connectivity = reshape (connectivity, sp.nsh_max, msh.nel);

clear conn_u conn_v conn_w connectivity

  ucp = sp.ndof_dir(1);
  vcp = sp.ndof_dir(2); 
  wcp = sp.ndof_dir(3);
  if (isfield (msh, 'boundary'))
    msh_bnd = msh.boundary;
    w_bnd{1} = sp.weights(1,:,:);
    w_bnd{2} = sp.weights(end,:,:);
    w_bnd{3} = sp.weights(:,1,:);
    w_bnd{4} = sp.weights(:,end,:);
    w_bnd{5} = sp.weights(:,:,1);
    w_bnd{6} = sp.weights(:,:,end);
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
        sp_bspline_2d_param (sp.knots(ind), sp.degree(ind), msh_bnd(iside));

      boundary = rmfield (bnd_iside, 'shape_function_gradients');
      boundary = bsp_2_nrb_2d__ (boundary, msh.boundary(iside), w_bnd{iside});

      boundary.dofs = sub2ind ([ucp, vcp, wcp], uidx, vidx, widx);

      sp.boundary(iside) = boundary;
    end 
  else
    sp.boundary = [];
  end

  sp = class (sp, 'sp_nurbs_3d');

end
