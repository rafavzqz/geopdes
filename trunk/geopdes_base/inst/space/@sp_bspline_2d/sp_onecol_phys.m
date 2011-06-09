% SP_EVALUATE_COL: compute the basis functions in one column of the mesh.
%
%     sp = sp_evaluate_col (space, msh, colnum, 'option1', value1, ...)
%
% INPUTS:
%     
%     sp:     class defining the space of discrete functions (see sp_bspline_2d)
%     msh:    msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization
%     colnum: number of the fixed element in the first parametric direction
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      true       |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%
%    FIELD_NAME      (SIZE)                      DESCRIPTION
%    ncomp           (scalar)                          number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                          total number of degrees of freedom
%    ndof_dir        (1 x 2 vector)                    degrees of freedom along each direction
%    nsh_max         (scalar)                          maximum number of shape functions per element
%    nsh             (1 x msh.nelv vector)             actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nelv vector)       indices of basis functions that do not vanish in each element
%    shape_functions (msh.nqn x nsh_max x msh.nelv)    basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%                        (2 x msh.nqn x nsh_max x msh.nelv) basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%                        (2 x 2 x msh.nqn x nsh_max x msh.nelv) basis function hessians evaluated at each quadrature node in each element
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

function [sp, elem_list] = sp_onecol_phys (space, msh, colnum, varargin)

value = true;
gradient = true;
hessian = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    else
      error ('sp_evaluate_col: unknown option %s', varargin {ii});
    end
  end
end

sp = sp_onecol_param (space, msh, colnum, varargin{:});

elem_list = colnum + msh.nelu*(0:msh.nelv-1);

if (gradient || hessian)
  if (gradient)
    JinvT = geopdes_invT__ (msh.geo_map_jac(:,:,:,elem_list));
    JinvT = reshape (JinvT, [2, 2, msh.nqn, msh.nelv]);
    sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
  end

  if (hessian && isfield (msh, 'geo_map_der2'))
    xu = squeeze (msh.geo_map_jac(1,1,:,elem_list));
    xv = squeeze (msh.geo_map_jac(1,2,:,elem_list)); 
    yu = squeeze (msh.geo_map_jac(2,1,:,elem_list)); 
    yv = squeeze (msh.geo_map_jac(2,2,:,elem_list)); 

    xuu = reshape (msh.geo_map_der2(1, 1, 1, :, elem_list), [], msh.nelv);
    yuu = reshape (msh.geo_map_der2(2, 1, 1, :, elem_list), [], msh.nelv);
    xuv = reshape (msh.geo_map_der2(1, 1, 2, :, elem_list), [], msh.nelv);
    yuv = reshape (msh.geo_map_der2(2, 2, 1, :, elem_list), [], msh.nelv);
    xvv = reshape (msh.geo_map_der2(1, 2, 2, :, elem_list), [], msh.nelv);
    yvv = reshape (msh.geo_map_der2(2, 2, 2, :, elem_list), [], msh.nelv);

    [uxx, uxy, uyy, vxx, vxy, vyy] = ...
            der2_inv_map__ (xu, xv, yu, yv, xuu, xuv, xvv, yuu, yuv, yvv);

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
    
    clear shape_function_hessians
    clear bu bv buu buv bvv xu xv yu yv xuu xuv xvv yuu yuv yvv ...
             uxx uxy uyy vxx vxy vyy bxx bxy byy
  end
  clear shg_u shg_v
end

clear shp_u shp_v

end
