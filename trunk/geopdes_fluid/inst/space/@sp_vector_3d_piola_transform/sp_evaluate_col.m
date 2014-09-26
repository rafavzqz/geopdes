% SP_EVALUATE_COL: compute the basis functions in one column of the mesh.
%
%     sp = sp_evaluate_col (space, msh_col, 'option1', value1, ...)
%
% INPUTS:
%     
%    sp:      object defining the space of discrete functions (see sp_vector_3d_piola_transform)
%    msh_col: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_3d/msh_evaluate_col)
%    option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            divergence |      false      |  compute shape_function_divs
%            curl       |      false      |  compute shape_function_curls
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                     DESCRIPTION
%    ncomp           (scalar)                                   number of components of the functions of the space (actually, 3)
%    ndof            (scalar)                                   total number of degrees of freedom
%    ndof_dir        (3 x 3 matrix)                             for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                                   maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)                   actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)             indices of basis functions that do not vanish in each element
%    shape_functions (3 x msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%               (3 x 3 x msh_col.nqn x nsh_max x msh_col.nel)   basis function gradients evaluated at each quadrature node in each element
%    shape_function_divs (msh_col.nqn x nsh_max x msh_col.nel)  basis function gradients evaluated at each quadrature node in each element
%    shape_function_curls (msh_col.nqn x nsh_max x msh_col.nel) basis function gradients evaluated at each quadrature node in each element
%
% Copyright (C) 2013 Elena Bulgarello, Carlo de Falco, Sara Frizziero
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

function sp = sp_evaluate_col (space, msh, varargin)

value = true;
gradient = false;
divergence = false;
curl = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'curl'))
      curl = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'divergence'))
      divergence = varargin {ii+1};
    else
      error ('sp_evaluate_col: unknown option %s', varargin {ii});
    end
  end
end

first_der = gradient || divergence || curl;
sp1_col = sp_evaluate_col_param (space.sp1, msh, 'value', true, 'gradient', first_der);
sp2_col = sp_evaluate_col_param (space.sp2, msh, 'value', true, 'gradient', first_der);
sp3_col = sp_evaluate_col_param (space.sp3, msh, 'value', true, 'gradient', first_der);


ndof     = sp1_col.ndof + sp2_col.ndof + sp3_col.ndof;
ndof_dir = [sp1_col.ndof_dir; sp2_col.ndof_dir; sp3_col.ndof_dir];
nsh      = sp1_col.nsh(:)' + sp2_col.nsh(:)' + sp3_col.nsh(:)';

connectivity = [sp1_col.connectivity; ...
                         sp2_col.connectivity+sp1_col.ndof; ...
                         sp3_col.connectivity + sp1_col.ndof + sp2_col.ndof];

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
            'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
            'ncomp', 3);

% We start computing the shape functions and gradients in the reference domain
shape_functions = zeros (3, msh.nqn, sp.nsh_max, msh.nel);
shape_functions(1,:,1:sp1_col.nsh_max,:)            = sp1_col.shape_functions;
shape_functions(2,:,sp1_col.nsh_max+(1:sp2_col.nsh_max),:) = sp2_col.shape_functions;
shape_functions(3,:,sp1_col.nsh_max+sp2_col.nsh_max+1:sp.nsh_max,:) = sp3_col.shape_functions;

if (gradient || curl)
  shape_function_gradients = zeros (3, 3, msh.nqn, sp.nsh_max, msh.nel);
  shape_function_gradients(1,:,:,1:sp1_col.nsh_max,:) = sp1_col.shape_function_gradients;
  shape_function_gradients(2,:,:, sp1_col.nsh_max+(1:sp2_col.nsh_max),:) = sp2_col.shape_function_gradients;
  shape_function_gradients(3,:,:, sp1_col.nsh_max+sp2_col.nsh_max+1:sp.nsh_max,:) = sp3_col.shape_function_gradients; 
end

% From here we apply the transformation
jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, msh.nel);

if (value)
  sp.shape_functions = geopdes_prod__ (msh.geo_map_jac, shape_functions);
  jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
  sp.shape_functions = bsxfun (@rdivide, sp.shape_functions, jacdet);
  jacdet = reshape (jacdet, msh.nqn, msh.nel);
end

if (divergence && ~(gradient || curl))
  sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);

  sp.shape_function_divs(:, 1:sp1_col.nsh_max, :) = reshape ...
     (sp1_col.shape_function_gradients(1,:,:,:), msh.nqn, sp1_col.nsh_max, msh.nel);
  sp.shape_function_divs(:, sp1_col.nsh_max+(1:sp2_col.nsh_max), :) = reshape ...
     (sp2_col.shape_function_gradients(2,:,:,:), msh.nqn, sp2_col.nsh_max, msh.nel);
  sp.shape_function_divs(:, sp1_col.nsh_max+sp2_col.nsh_max+1:sp.nsh_max, :) = reshape ...
     (sp3_col.shape_function_gradients(3,:,:,:), msh.nqn, sp3_col.nsh_max, msh.nel);

  jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
  sp.shape_function_divs = bsxfun (@rdivide, sp.shape_function_divs, jacdet);
  jacdet = reshape (jacdet, msh.nqn, msh.nel);
end

if (gradient || curl)
  shape_fun_grads = zeros (2, 2, msh.nqn, sp.nsh_max, msh.nel);

  xu = reshape (msh.geo_map_jac(1, 1, :, :), msh.nqn, msh.nel);
  xv = reshape (msh.geo_map_jac(1, 2, :, :), msh.nqn, msh.nel);
  xw = reshape (msh.geo_map_jac(1, 3, :, :), msh.nqn, msh.nel);
  yu = reshape (msh.geo_map_jac(2, 1, :, :), msh.nqn, msh.nel);
  yv = reshape (msh.geo_map_jac(2, 2, :, :), msh.nqn, msh.nel);
  yw = reshape (msh.geo_map_jac(2, 3, :, :), msh.nqn, msh.nel);
  zu = reshape (msh.geo_map_jac(3, 1, :, :), msh.nqn, msh.nel);
  zv = reshape (msh.geo_map_jac(3, 2, :, :), msh.nqn, msh.nel);
  zw = reshape (msh.geo_map_jac(3, 3, :, :), msh.nqn, msh.nel);

  xuu = reshape (msh.geo_map_der2(1, 1, 1, :, :), [], msh.nel);
  yuu = reshape (msh.geo_map_der2(2, 1, 1, :, :), [], msh.nel);
  zuu = reshape (msh.geo_map_der2(3, 1, 1, :, :), [], msh.nel);
  xuv = reshape (msh.geo_map_der2(1, 1, 2, :, :), [], msh.nel);
  yuv = reshape (msh.geo_map_der2(2, 2, 1, :, :), [], msh.nel);
  zuv = reshape (msh.geo_map_der2(3, 2, 1, :, :), [], msh.nel);
  xvv = reshape (msh.geo_map_der2(1, 2, 2, :, :), [], msh.nel);
  yvv = reshape (msh.geo_map_der2(2, 2, 2, :, :), [], msh.nel);
  zvv = reshape (msh.geo_map_der2(3, 2, 2, :, :), [], msh.nel);
  xuw= reshape (msh.geo_map_der2(1, 1, 3, :, :), [], msh.nel);
  yuw= reshape (msh.geo_map_der2(2, 1, 3, :, :), [], msh.nel);
  zuw= reshape (msh.geo_map_der2(3, 3, 1, :, :), [], msh.nel);
  xvw= reshape (msh.geo_map_der2(1, 2, 3, :, :), [], msh.nel);
  yvw= reshape (msh.geo_map_der2(2, 2, 3, :, :), [], msh.nel);
  zvw= reshape (msh.geo_map_der2(3, 3, 2, :, :), [], msh.nel);
  xww= reshape (msh.geo_map_der2(1, 3, 3, :, :), [], msh.nel);
  yww= reshape (msh.geo_map_der2(2, 3, 3, :, :), [], msh.nel);
  zww= reshape (msh.geo_map_der2(3, 3, 3, :, :), [], msh.nel);
 
  for ii=1:sp.nsh_max
    wh  = reshape (shape_functions(1, :, ii, :), msh.nqn, msh.nel);
    zh  = reshape (shape_functions(2, :, ii, :), msh.nqn, msh.nel);
    lh  = reshape (shape_functions(3, :, ii, :), msh.nqn, msh.nel);
    whu = reshape (shape_function_gradients(1, 1, :, ii, :), msh.nqn, msh.nel);
    whv = reshape (shape_function_gradients(1, 2, :, ii, :), msh.nqn, msh.nel);
    whw = reshape (shape_function_gradients(1, 3, :, ii, :), msh.nqn, msh.nel);
    zhu = reshape (shape_function_gradients(2, 1, :, ii, :), msh.nqn, msh.nel);
    zhv = reshape (shape_function_gradients(2, 2, :, ii, :), msh.nqn, msh.nel);
    zhw = reshape (shape_function_gradients(2, 3, :, ii, :), msh.nqn, msh.nel);
    lhu = reshape (shape_function_gradients(3, 1, :, ii, :), msh.nqn, msh.nel);
    lhv = reshape (shape_function_gradients(3, 2, :, ii, :), msh.nqn, msh.nel);
    lhw = reshape (shape_function_gradients(3, 3, :, ii, :), msh.nqn, msh.nel);
    [wx, wy, wz, zx, zy, zz, lx, ly, lz] = piola_transform_grad_3d__ (...
                        xuu, yuu, zuu, xuv, yuv, zuv, xuw, yuw, zuw, xvv, yvv, zvv, ...
     xvw, yvw, zvw, xww, yww, zww, xu, xv, xw, yu, yv, yw, zu, zv, zw,...
                        wh, zh, lh, whu, whv, whw, zhu, zhv, zhw, lhu, lhv, lhw);
                        
    shape_fun_grads(1, 1, :, ii, :) = wx;
    shape_fun_grads(1, 2, :, ii, :) = wy;
    shape_fun_grads(1, 3, :, ii, :) = wz;
    shape_fun_grads(2, 1, :, ii, :) = zx;
    shape_fun_grads(2, 2, :, ii, :) = zy;
    shape_fun_grads(2, 3, :, ii, :) = zz;
    shape_fun_grads(3, 1, :, ii, :) = lx;
    shape_fun_grads(3, 2, :, ii, :) = ly;
    shape_fun_grads(3, 3, :, ii, :) = lz;
  end

  if (gradient)
    sp.shape_function_gradients = shape_fun_grads;
  end

  if (divergence)
    sp.shape_function_divs = reshape (shape_fun_grads(1,1,:,:,:) + ...
				      shape_fun_grads(2,2,:,:,:) + shape_fun_grads(3,3,:,:,:), ...
                                      msh.nqn, sp.nsh_max, msh.nel);
  end

  if (curl)
    shape_fun_curls = zeros (3, msh.nqn, sp.nsh_max, msh.nel);
    shape_fun_curls(1,:,:,:) = ...
      reshape (shape_fun_grads(3,2,:,:,:) - shape_fun_grads(2,3,:,:,:), ...
                                       1, msh.nqn, sp.nsh_max, msh.nel);
    shape_fun_curls(2,:,:,:) = ...
      reshape (shape_fun_grads(1,3,:,:,:) - shape_fun_grads(3,1,:,:,:), ...
                                       1, msh.nqn, sp.nsh_max, msh.nel);
    shape_fun_curls(3,:,:,:) = ...
      reshape (shape_fun_grads(2,1,:,:,:) - shape_fun_grads(1,2,:,:,:), ...
                                       1, msh.nqn, sp.nsh_max, msh.nel);
    sp.shape_function_curls = shape_fun_curls;                                   
                                       
                                       
  end
end

end
