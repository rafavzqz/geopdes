% SP_BSPLINE_2D_PHYS: Construct a tensor-product space of B-Splines on the physical domain in 2D.
%
%     sp = sp_bspline_2d_phys (knots, degree, msh, 'option1', value1, ...)
%
% INPUTS:
%     
%     knots:  open knot vector    
%     degree: b-spline polynomial degree
%     msh:    structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+-----------------------------------
%            gradient   |      true       |  compute shape_function_gradients
%            hessians   |      false      |  compute shape_function_hessians
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
%        shape_function_hessians
%                        (2 x 2 x msh.nqn x nsh_max x msh.nel) basis function hessians evaluated at each quadrature node in each element
%        boundary        (1 x 4 struct array)              struct array representing the space of traces of basis functions on each edge
%        spfun           (function handle)                 function to evaluate an element of the discrete function space, given the Fourier 
%                                                          coefficients and a set of points in the parametric space
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

function sp = sp_bspline_2d_phys (knots, degree, msh, varargin)

sp = sp_bspline_2d_param (knots, degree, msh, varargin{:});


if (isfield (sp, 'shape_function_hessians'))
  if (isfield (msh, 'geo_map_der2'))
    
    xu = squeeze (msh.geo_map_jac(1,1,:,:));
    xv = squeeze (msh.geo_map_jac(1,2,:,:)); 
    yu = squeeze (msh.geo_map_jac(2,1,:,:)); 
    yv = squeeze (msh.geo_map_jac(2,2,:,:)); 
    
    xuu = reshape (msh.geo_map_der2(1, 1, 1, :, :), [], msh.nel);
    yuu = reshape (msh.geo_map_der2(2, 1, 1, :, :), [], msh.nel);
    xuv = reshape (msh.geo_map_der2(1, 1, 2, :, :), [], msh.nel);
    yuv = reshape (msh.geo_map_der2(2, 2, 1, :, :), [], msh.nel);
    xvv = reshape (msh.geo_map_der2(1, 2, 2, :, :), [], msh.nel);
    yvv = reshape (msh.geo_map_der2(2, 2, 2, :, :), [], msh.nel);
    
    [uxx, uxy, uyy, vxx, vxy, vyy] = der2_inv_map__ (xu, xv, yu, yv, xuu, xuv, xvv, yuu, yuv, yvv);

    for ii=1:sp.nsh_max   
      bu = squeeze (sp.shape_function_gradients(1,:,ii,:));
      bv = squeeze (sp.shape_function_gradients(2,:,ii,:));
      buu = squeeze (sp.shape_function_hessians(1,1,:,ii,:));
      buv = squeeze (sp.shape_function_hessians(1,2,:,ii,:));
      bvv = squeeze (sp.shape_function_hessians(2,2,:,ii,:));
      
      [bxx, bxy, byy] = der2_basisfun_phys__ (xu(:), xv(:), yu(:), yv(:), uxx(:), uxy(:), uyy(:), vxx(:), vxy(:), vyy(:), buu(:), buv(:), bvv(:), bu(:), bv(:));
      sh = size (sp.shape_function_hessians(1,1,:,ii,:));
      shape_function_hessians(1,1,:,ii,:) = reshape (bxx, sh);
      shape_function_hessians(1,2,:,ii,:) = reshape (bxy, sh);
      shape_function_hessians(2,1,:,ii,:) = reshape (bxy, sh);
      shape_function_hessians(2,2,:,ii,:) = reshape (byy, sh);
    end  
    sp.shape_function_hessians = shape_function_hessians;
  else
    rmfield (sp, 'shape_function_hessians');
  end
end

if (isfield (sp, 'shape_function_gradients'))
  JinvT = geopdes_invT__ (msh.geo_map_jac);
  JinvT = reshape (JinvT, [2, 2, msh.nqn, msh.nel]);
  shape_fun_grads = reshape (sp.shape_function_gradients, [2, msh.nqn, sp.nsh_max, msh.nel]);
  sp.shape_function_gradients = geopdes_prod__ (JinvT, shape_fun_grads);
end

sp.spfun  = @(MSH) sp_bspline_2d_phys (knots, degree, MSH, varargin{:});

end

