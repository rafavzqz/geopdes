% SP_BSPLINE_3D_PARAM: Construct a tensor-product space of B-Splines on the parametric domain in 3D.
%                      This function is not usually meant to be invoked directly by the user but rather
%                      through sp_bspline_3d_phys.

%     sp = sp_bspline_3d_param (knots, degree, msh, 'option1', value1, ...)
%
% INPUTS:
%     
%     knots:   open knot vector    
%     degree:  b-spline polynomial degree
%     msh:     msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization  
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
%        ndof_dir        (1 x 3 vector)                    degrees of freedom along each direction (only for tensor product spaces)
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)              actual number of shape functions per each element
%        connectivity    (nsh_max x msh.nel vector)        indices of basis functions that do not vanish in each element
%        shape_functions (msh.nqn x nsh_max x msh.nel)     basis functions evaluated at each quadrature node in each element
%        shape_function_gradients
%                        (3 x msh.nqn x nsh_max x msh.nel) basis function gradients evaluated at each quadrature node in each element
%        boundary        (1 x 6 struct array)              struct array representing the space of traces of basis functions on each face
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

function sp = sp_bspline_3d_param (knots, degree, msh, varargin)

  gradient = true; 
  if (~isempty (varargin))
    if (~rem (length (varargin), 2) == 0)
      error ('sp_bspline_3d_param: options must be passed in the [option, value] format');
    end
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin {ii}, 'gradient'))
        gradient = varargin {ii+1};
      else
        error ('sp_bspline_3d_param: unknown option %s', varargin {ii});
      end
    end
  end


  nodes = msh.qn;
  spu = sp_bspline_1d_param (knots{1}, degree(1), nodes{1});
  spv = sp_bspline_1d_param (knots{2}, degree(2), nodes{2});
  spw = sp_bspline_1d_param (knots{3}, degree(3), nodes{3});

  nsh_max  = spu.nsh_max * spv.nsh_max * spw.nsh_max;
  [NSHV, NSHU, NSHW] = meshgrid (spv.nsh, spu.nsh, spw.nsh);
  nsh  = reshape (NSHU .* NSHV .* NSHW, 1, []);
  ndof = spu.ndof * spv.ndof * spw.ndof;
  ndof_dir = [spu.ndof, spv.ndof, spw.ndof];

  nelu = size (nodes {1}, 2);
  nelv = size (nodes {2}, 2);
  nelw = size (nodes {3}, 2);
  nel  = nelu * nelv * nelw;

  nqnu = size (nodes {1}, 1);
  nqnv = size (nodes {2}, 1);
  nqnw = size (nodes {3}, 1);
  nqn  = nqnu * nqnv * nqnw;

  conn_u = reshape (spu.connectivity, spu.nsh_max, 1, 1, nelu, 1, 1);
  conn_u = repmat  (conn_u, [1, spv.nsh_max, spw.nsh_max, 1, nelv, nelw]);
  conn_u = reshape (conn_u, [], nel);

  conn_v = reshape (spv.connectivity, 1, spv.nsh_max, 1, 1, nelv, 1);
  conn_v = repmat  (conn_v, [spu.nsh_max, 1, spw.nsh_max, nelu, 1, nelw]);
  conn_v = reshape (conn_v, [], nel);

  conn_w = reshape (spw.connectivity, 1, 1, spw.nsh_max, 1, 1, nelw);
  conn_w = repmat  (conn_w, [spu.nsh_max, spv.nsh_max, 1, nelu, nelv, 1]);
  conn_w = reshape (conn_w, [], nel);

  connectivity = zeros (nsh_max, nel);
  indices = conn_u ~= 0 & conn_v ~= 0 & conn_w~=0;
  connectivity(indices) = ...
     sub2ind ([spu.ndof, spv.ndof, spw.ndof], conn_u(indices), conn_v(indices), conn_w(indices));
  connectivity = reshape (connectivity, nsh_max, nel);

  clear conn_u conn_v conn_w

  shp_u = reshape (spu.shape_functions, nqnu, 1, 1, spu.nsh_max, 1, 1, nelu, 1, 1);
  shp_u = repmat  (shp_u, [1, nqnv, nqnw, 1, spv.nsh_max, spw.nsh_max, 1, nelv, nelw]);
  shp_u = reshape (shp_u, nqn, nsh_max, nel);

  shp_v = reshape (spv.shape_functions, 1, nqnv, 1, 1, spv.nsh_max, 1, 1, nelv, 1);
  shp_v = repmat  (shp_v, [nqnu, 1, nqnw, spu.nsh_max, 1, spw.nsh_max, nelu, 1, nelw]);
  shp_v = reshape (shp_v, nqn, nsh_max, nel);

  shp_w = reshape (spw.shape_functions, 1, 1, nqnw, 1, 1, spw.nsh_max, 1, 1, nelw);
  shp_w = repmat  (shp_w, [nqnu, nqnv, 1, spu.nsh_max, spv.nsh_max, 1, nelu, nelv, 1]);
  shp_w = reshape (shp_w, nqn, nsh_max, nel);

  shape_functions = shp_u .* shp_v .* shp_w;

  sp = struct('nsh_max', nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
	      'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
	      'shape_functions', shape_functions, ...
	      'ncomp', 1);

  if (gradient)
    shg_u = reshape (spu.shape_function_gradients, nqnu, 1, 1, spu.nsh_max, 1, 1, nelu, 1, 1);
    shg_u = repmat  (shg_u, [1, nqnv, nqnw, 1, spv.nsh_max, spw.nsh_max, 1, nelv, nelw]);
    shg_u = reshape (shg_u, nqn, nsh_max, nel);
    
    shg_v = reshape (spv.shape_function_gradients, 1, nqnv, 1, 1, spv.nsh_max, 1, 1, nelv, 1);
    shg_v = repmat  (shg_v, [nqnu, 1, nqnw, spu.nsh_max, 1, spw.nsh_max, nelu, 1, nelw]);
    shg_v = reshape (shg_v, nqn, nsh_max, nel);
    
    shg_w = reshape (spw.shape_function_gradients, 1, 1, nqnw, 1, 1, spw.nsh_max, 1, 1, nelw);
    shg_w = repmat  (shg_w, [nqnu, nqnv, 1, spu.nsh_max, spv.nsh_max, 1, nelu, nelv, 1]);
    shg_w = reshape (shg_w, nqn, nsh_max, nel);
    
    shape_function_gradients(1,:,:,:) = shg_u .* shp_v .* shp_w ;
    shape_function_gradients(2,:,:,:) = shp_u .* shg_v .* shp_w ;
    shape_function_gradients(3,:,:,:) = shp_u .* shp_v .* shg_w ;

    clear  shg_u shg_v shg_w
    sp.shape_function_gradients = shape_function_gradients;
  end  
  
  clear shp_u shp_v shp_w
    
  ucp = ndof_dir(1);
  vcp = ndof_dir(2);
  wcp = ndof_dir(3);
  
  
  bpoints = msh.boundary;
  for iside = msh.boundary_list
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
    sp_bspline_2d_param (knots(ind), degree(ind), bpoints(iside));

    boundary      = rmfield (bnd_iside, 'shape_function_gradients');
    boundary.dofs = sub2ind ([ucp, vcp, wcp], uidx, vidx, widx);

    sp.boundary(iside) = boundary;
  end 
  
  sp.spfun = @(MSH) sp_bspline_3d_param (knots, degree, MSH);

end

%!test
%! knots = {[0 0 0 .5 1 1 1], [0 0 0 .25 .75 1 1 1], [0 0 0 .25 .75 1 1 1]};
%! degree = [2 2 2];
%! nodes{1} = [0 .1 .2; .6 .7 .8]';
%! nodes{2} = [.1 .2; .3 .7; .8 .9]';
%! nodes{3} = [.1 .2; .3 .7; .8 .9]';
%! msh.qn = nodes;
%! sp = sp_bspline_3d_param (knots, degree, msh);
%! assert (norm (sum (sp.shape_functions, 2)(:) - 1, inf), 0, 1e-10)
%! assert (norm (sum (sp.shape_function_gradients, 3)(:), inf), 0, 1e-10)
