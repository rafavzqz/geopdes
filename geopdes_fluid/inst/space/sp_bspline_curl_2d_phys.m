% SP_BSPLINE_CURL_2D_PHYS:
%
%     sp = sp_bspline_curl_2d_phys (knots, degree, msh)
%
% INPUTS:
%     
%     knots:  open knot vector    
%     degree: b-spline polynomial degree
%     msh:    structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
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
%        boundary        (1 x 4 struct array)              *** INCOMPLETE IN CURRENT IMPLEMENTATION ***
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

function sp = sp_bspline_curl_2d_phys (knots, degree, msh)

spb = sp_bspline_2d_phys (knots, degree, msh, 'hessian', true);

sp.ncomp        = 2;
sp.ndof         = spb.ndof;
sp.nsh_max      = spb.nsh_max;
sp.nsh          = spb.nsh;
sp.connectivity = spb.connectivity;

sp.shape_functions(1,:,:,:) =  squeeze (spb.shape_function_gradients (2,:,:,:));
sp.shape_functions(2,:,:,:) = -squeeze (spb.shape_function_gradients (1,:,:,:));

sp.shape_function_gradients(1,1,:,:,:) =  spb.shape_function_hessians (2,1,:,:,:);
sp.shape_function_gradients(1,2,:,:,:) =  spb.shape_function_hessians (2,2,:,:,:);
sp.shape_function_gradients(2,1,:,:,:) = -spb.shape_function_hessians (1,1,:,:,:);
sp.shape_function_gradients(2,2,:,:,:) = -spb.shape_function_hessians (1,2,:,:,:);

if (isfield (msh, 'boundary'))
  mcp = spb.ndof_dir(1);
  ncp = spb.ndof_dir(2);

  sp.boundary(1).dofs = union (sub2ind ([mcp, ncp], ones(1,ncp), 1:ncp), sub2ind ([mcp, ncp], 2*ones(1,ncp), 1:ncp));
  sp.boundary(2).dofs = union (sub2ind ([mcp, ncp], mcp*ones(1,ncp), 1:ncp), sub2ind ([mcp, ncp], (mcp-1)*ones(1,ncp), 1:ncp));
  sp.boundary(3).dofs = union (sub2ind ([mcp, ncp], 1:mcp, ones(1,mcp)), sub2ind ([mcp, ncp], 1:mcp, 2*ones(1,mcp)));
  sp.boundary(4).dofs = union (sub2ind ([mcp, ncp], 1:mcp, ncp*ones(1,mcp)), sub2ind ([mcp, ncp], 1:mcp, (ncp-1)*ones(1,mcp)));

%  sp.boundary(1).dofs = spb.boundary(1).dofs; 
%  sp.boundary(2).dofs = spb.boundary(2).dofs; 
%  sp.boundary(3).dofs = spb.boundary(3).dofs; 
%  sp.boundary(4).dofs = spb.boundary(4).dofs; 
end

sp.spfun  = @(MSH) sp_bspline_curl_2d_phys (knots, degree, MSH);

end

