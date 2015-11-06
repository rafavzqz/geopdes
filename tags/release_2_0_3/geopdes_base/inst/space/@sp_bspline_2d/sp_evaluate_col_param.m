% SP_EVALUATE_COL_PARAM: compute the basis functions, in the parametric domain, in one column of the mesh.
%
%     sp = sp_evaluate_col_param (space, msh_col, 'option1', value1, ...)
%
% INPUTS:
%
%    space:   object defining the space of discrete functions (see sp_bspline_2d)
%    msh_col: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_2d/msh_evaluate_col)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                 DESCRIPTION
%    ncomp           (scalar)                               number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                               total number of degrees of freedom
%    ndof_dir        (1 x 2 vector)                         degrees of freedom along each direction
%    nsh_max         (scalar)                               maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%            (2 x msh_col.nqn x nsh_max x msh_col.nel)      basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%            (2 x 2 x msh_col.nqn x nsh_max x msh_col.nel)  basis function hessians evaluated at each quadrature node in each element
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

function sp = sp_evaluate_col_param (space, msh, varargin)

value = true;
gradient = false;
hessian = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col_param: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    else
      error ('sp_evaluate_col_param: unknown option %s', varargin {ii});
    end
  end
end

spu = space.spu;
spv = space.spv;

ndof = spu.ndof * spv.ndof;
ndof_dir = [spu.ndof, spv.ndof];

nsh  = spu.nsh(msh.colnum) * spv.nsh;
nsh  = nsh(:)';

conn_u = reshape (spu.connectivity(:,msh.colnum), spu.nsh_max, 1, 1);
conn_u = repmat  (conn_u, [1, spv.nsh_max, msh.nel]);
conn_u = reshape (conn_u, [], msh.nel);

conn_v = reshape (spv.connectivity, 1, spv.nsh_max, msh.nel_dir(2));
conn_v = repmat  (conn_v, [spu.nsh_max, 1, 1]);
conn_v = reshape (conn_v, [], msh.nel);

connectivity = zeros (space.nsh_max, msh.nel);
indices = (conn_u ~= 0) & (conn_v ~= 0);
connectivity(indices) = ...
  sub2ind ([spu.ndof, spv.ndof], conn_u(indices), conn_v(indices));
connectivity = reshape (connectivity, space.nsh_max, msh.nel);

clear conn_u conn_v

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
            'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
            'ncomp', 1);

shp_u = reshape (spu.shape_functions(:, :, msh.colnum), ...
               msh.nqn_dir(1), 1, spu.nsh_max, 1, 1);  %% one column only
shp_u = repmat  (shp_u, [1, msh.nqn_dir(2), 1, spv.nsh_max, msh.nel]);
shp_u = reshape (shp_u, msh.nqn, space.nsh_max, msh.nel);

shp_v = reshape (spv.shape_functions, 1, msh.nqn_dir(2), 1, spv.nsh_max, msh.nel);
shp_v = repmat  (shp_v, [msh.nqn_dir(1), 1, spu.nsh_max, 1, 1]);
shp_v = reshape (shp_v, msh.nqn, space.nsh_max, msh.nel);

if (value)
  sp.shape_functions = shp_u .* shp_v ;
end

if (gradient || hessian)
  shg_u = reshape (spu.shape_function_gradients(:,:,msh.colnum), ...
                   msh.nqn_dir(1), 1, spu.nsh_max, 1, 1);  %% one column only
  shg_u = repmat  (shg_u, [1, msh.nqn_dir(2), 1, spv.nsh_max, msh.nel]);
  shg_u = reshape (shg_u, msh.nqn, space.nsh_max, msh.nel);

  shg_v = reshape (spv.shape_function_gradients, ...
                   1, msh.nqn_dir(2), 1, spv.nsh_max, msh.nel);
  shg_v = repmat  (shg_v, [msh.nqn_dir(1), 1, spu.nsh_max, 1, 1]);
  shg_v = reshape (shg_v, msh.nqn, space.nsh_max, msh.nel);

  sp.shape_function_gradients(1,:,:,:) = shg_u .* shp_v ;
  sp.shape_function_gradients(2,:,:,:) = shp_u .* shg_v ;

  sp.shape_function_gradients = reshape (sp.shape_function_gradients, ...
                                2, msh.nqn, sp.nsh_max, msh.nel);

  if (hessian && isfield (msh, 'geo_map_der2'))
    shh_uu = reshape (spu.shape_function_hessians(:,:,msh.colnum), ...
                      msh.nqn_dir(1), 1, spu.nsh_max, 1, 1);
    shh_uu = repmat  (shh_uu, [1, msh.nqn_dir(2), 1, spv.nsh_max, msh.nel]);
    shh_uu = reshape (shh_uu, msh.nqn, sp.nsh_max, msh.nel);

    shh_vv = reshape (spv.shape_function_hessians, 1, msh.nqn_dir(2), 1, spv.nsh_max, msh.nel);
    shh_vv = repmat  (shh_vv, [msh.nqn_dir(1), 1, spu.nsh_max, 1, 1]);
    shh_vv = reshape (shh_vv, msh.nqn, sp.nsh_max, msh.nel);

    shape_function_hessians(1,1,:,:,:) = shh_uu .* shp_v ;
    shape_function_hessians(1,2,:,:,:) = shg_u  .* shg_v ;
    shape_function_hessians(2,1,:,:,:) = shape_function_hessians(1,2,:,:,:);
    shape_function_hessians(2,2,:,:,:) = shp_u .* shh_vv ;
    sp.shape_function_hessians = shape_function_hessians;

    clear shh_uu shh_vv shape_function_hessians
  end
  clear shg_u shg_v
end

clear shp_u shp_v

end
