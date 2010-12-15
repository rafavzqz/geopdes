% SP_BSPLINE_2D_PARAM: Construct a tensor-product space of B-Splines on the parametric domain in 2D.
%                      This function is not usually meant to be invoked directly by the user but rather
%                      through sp_bspline_2d_phys.
%
%     sp = sp_bspline_2d_param (knots, degree, msh, 'option1', value1, ...)
%
% INPUTS:
%     
%     knots:     open knot vector    
%     degree:    b-spline polynomial degree
%     msh:       msh structure containing (in the field msh.qn) the points 
%                along each parametric direction in the parametric 
%                domain at which to evaluate, i.e. quadrature points 
%                or points for visualization  
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+-----------
%            gradient   |      true       |  compute shape_function_gradients
%
% OUTPUT:
%
%    sp: structure representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        ndof            (scalar)                          total number of degrees of freedom
%        ndof_dir        (1 x 2 vector)                    degrees of freedom along each direction (only for tensor product spaces)
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)              actual number of shape functions per each element
%        connectivity    (nsh_max x msh.nel vector)        indices of basis functions that do not vanish in each element
%        shape_functions (msh.nqn x nsh_max x msh.nel)     basis functions evaluated at each quadrature node in each element
%        shape_function_gradients
%                        (2 x msh.nqn x nsh_max x msh.nel) basis function gradients evaluated at each quadrature node in each element
%        boundary        (1 x 4 struct array)              struct array representing the space of traces of basis functions on each edge
%
%   For more details, see the documentation
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function sp = sp_bspline_2d_param (knots, degree, msh, varargin)

  gradient = true; 
  if (~isempty (varargin))
    if (~rem (length (varargin), 2) == 0)
      error ('sp_bspline_2d_param: options must be passed in the [option, value] format');
    end
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin {ii}, 'gradient'))
        gradient = varargin {ii+1};
      else
        error ('sp_bspline_2d_param: unknown option %s', varargin {ii});
      end
    end
  end

  nodes = msh.qn;
  spu = sp_bspline_1d_param (knots{1}, degree(1), nodes{1});
  spv = sp_bspline_1d_param (knots{2}, degree(2), nodes{2});

  nsh_max  = spu.nsh_max * spv.nsh_max;
  nsh  = spu.nsh' * spv.nsh;
  nsh  = nsh(:)';
  ndof = spu.ndof * spv.ndof;
  ndof_dir = [spu.ndof, spv.ndof];

  nelu = size (nodes{1}, 2);
  nelv = size (nodes{2}, 2);
  nel  = nelu * nelv;

  nqnu = size (nodes{1}, 1);
  nqnv = size (nodes{2}, 1);
  nqn  = nqnu * nqnv;

  conn_u = reshape (spu.connectivity, spu.nsh_max, 1, nelu, 1);
  conn_u = repmat  (conn_u, [1, spv.nsh_max, 1, nelv]);

  conn_v = reshape (spv.connectivity, 1, spv.nsh_max, 1, nelv);
  conn_v = repmat  (conn_v, [spu.nsh_max, 1, nelu, 1]);

  indices = (conn_u ~= 0) & (conn_v ~= 0);
  connectivity = zeros (nsh_max, nel);
  connectivity(indices) = sub2ind ([spu.ndof, spv.ndof], conn_u(indices),conn_v(indices));
  connectivity = reshape (connectivity, [nsh_max, nel]);
  
  % magick trick to sort zeros to the end of each column as
  % suggested by b. abbot: http://octave.1599824.n4.nabble.com/Yet-Another-Vectorization-Quiz-td3088102.html#a3088146
  [ignore, indices] = sort (connectivity==0);
  connectivity = connectivity (ones (size (connectivity, 1), 1) * (0:size (connectivity, 2) - 1) * size (connectivity, 1) + indices);

  clear conn_u conn_v

  shp_u = reshape (spu.shape_functions, nqnu, 1, spu.nsh_max, 1, nelu, 1);
  shp_u = repmat  (shp_u, [1, nqnv, 1, spv.nsh_max, 1, nelv]);
  shp_u = reshape (shp_u, nqn, nsh_max, nel);

  shp_v = reshape (spv.shape_functions, 1, nqnv, 1, spv.nsh_max, 1, nelv);
  shp_v = repmat  (shp_v, [nqnu, 1, spu.nsh_max, 1, nelu, 1]);
  shp_v = reshape (shp_v, nqn, nsh_max, nel);

  shape_functions = shp_u .* shp_v ;

  sp = struct('nsh_max', nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
	      'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
	      'shape_functions', shape_functions, ...
	      'ncomp', 1);

  if (gradient)
  
    shg_u = reshape (spu.shape_function_gradients, nqnu, 1, spu.nsh_max, 1, nelu, 1);
    shg_u = repmat  (shg_u, [1, nqnv, 1, spv.nsh_max, 1, nelv]);
    shg_u = reshape (shg_u, nqn, nsh_max, nel);
    
    shg_v = reshape (spv.shape_function_gradients, 1, nqnv, 1, spv.nsh_max, 1, nelv);
    shg_v = repmat  (shg_v, [nqnu, 1, spu.nsh_max, 1, nelu, 1]);
    shg_v = reshape (shg_v, nqn, nsh_max, nel);
    
    shape_function_gradients(1,:,:,:) = shg_u .* shp_v ;
    shape_function_gradients(2,:,:,:) = shp_u .* shg_v ;
    
    clear shg_u shg_v
    sp.shape_function_gradients = shape_function_gradients;
    
  end
  
  clear shp_u shp_v
  
  mcp = ndof_dir(1);
  ncp = ndof_dir(2);

  if (isfield (msh, 'boundary'))
    for iside = 1:numel(msh.boundary)
      ind = mod (floor ((iside+1)/2), 2) + 1;
      bnodes = reshape (squeeze (msh.boundary(iside).quad_nodes(ind,:,:)), ...
                msh.boundary(iside).nqn, []);
      bnd_iside = sp_bspline_1d_param (knots{ind}, degree(ind), bnodes);
      boundary(iside) = rmfield (bnd_iside, 'shape_function_gradients');
    end

    boundary(1).dofs = sub2ind ([mcp, ncp], ones(1,ncp), 1:ncp);
    boundary(2).dofs = sub2ind ([mcp, ncp], mcp*ones(1,ncp), 1:ncp);
    boundary(3).dofs = sub2ind ([mcp, ncp], 1:mcp, ones(1,mcp));
    boundary(4).dofs = sub2ind ([mcp, ncp], 1:mcp, ncp*ones(1,mcp));

    sp.boundary = boundary;
  end
  sp.spfun = @(MSH) sp_bspline_2d_param (knots, degree, MSH);

end

%!test
%! knots = {[0 0 0 .5 1 1 1], [0 0 0 .25 .75 1 1 1]};
%! degree = [2 2];
%! nodes{1} = [0 .1 .2; .6 .7 .8]';
%! nodes{2} = [.1 .2; .3 .7; .8 .9]';
%! msh.qn = nodes; 
%! sp = sp_bspline_2d_param (knots, degree, msh);
%! assert (norm (sum (sp.shape_functions, 2)(:) - 1, inf), 0, 1e-10)
%! assert (norm (sum (sp.shape_function_gradients, 3)(:), inf), 0, 1e-10)
