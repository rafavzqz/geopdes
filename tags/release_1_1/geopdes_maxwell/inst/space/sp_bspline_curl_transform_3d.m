% SP_BSPLINE_CURL_TRANSFORM_3D: Construct a space of B-Splines mapping to the physical domain with a curl conserving transform.
%
%     sp = sp_bspline_curl_transform_3d (knots_u1, knots_u2, knots_u3, ...
%                     degree1, degree2, degree3, msh)
%
% INPUTS:
%
%     knots_u1: open knot vector for the first component of the vector
%     knots_u2: open knot vector for the second component of the vector
%     knots_u3: open knot vector for the third component of the vector
%     degree1:  b-spline polynomial degree for space of the first component.
%     degree2:  b-spline polynomial degree for space of the second component.
%     degree3:  b-spline polynomial degree for space of the third component.
%                The three spaces should be of the form 
%                     S^{p,q+1,r+1} x S^{p+1,q,r+1} x S^{p+1,q+1,r}
%     msh:      structure containing the domain partition and the quadrature rule (see msh_push_forward_3d)
%
% OUTPUT:
%
%    sp: structure representing the discrete function space, with the following fields:
%
%        FIELD_NAME      (SIZE)                            DESCRIPTION
%        ndof            (scalar)                          total number of degrees of freedom
%        nsh_max         (scalar)                          maximum number of shape functions per element
%        nsh             (1 x msh.nel vector)              actual number of shape functions per each element
%        connectivity    (nsh_max x msh.nel vector)        indices of basis functions that do not vanish in each element
%        shape_functions (3 x msh.nqn x nsh_max x msh.nel) vector-valued basis functions evaluated at each quadrature node in each element
%        shape_function_curls
%                        (3 x msh.nqn x nsh_max x msh.nel) basis functions curls evaluated at each quadrature node in each element
%        boundary        (1 x 6 struct array)              struct array representing the space of tangential traces of basis functions on each edge
%        spfun           (function handle)                 function to evaluate an element of the discrete function space, given the Fourier coefficients and a set of points in the parametric space
%
%   For more details, see the documentation
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010 Rafael Vazquez
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
  
function sp = sp_bspline_curl_transform_3d (knots_u1, knots_u2, knots_u3, ...
                         degree1, degree2, degree3, msh)


  sp1 = sp_bspline_3d_param (knots_u1, degree1, msh);
  sp2 = sp_bspline_3d_param (knots_u2, degree2, msh);
  sp3 = sp_bspline_3d_param (knots_u3, degree3, msh);
  
  sp  = sp_scalar_to_vector_3d (sp1, sp2, sp3, msh, 'curl', true, 'gradient', false);

% Compute the boundary field, which only contains tangential components
  if (isfield (msh, 'boundary'))
    sp.boundary = rmfield (sp.boundary, 'spfun');
    for iface = 1:6
      switch (iface)
       case {1,2}
        boundary.ncomp    = 3;
        boundary.ndof     = ...
                 sp2.boundary(iface).ndof + sp3.boundary(iface).ndof;
        boundary.ndof_dir = ...
                [sp2.boundary(iface).ndof_dir; sp3.boundary(iface).ndof_dir];
        boundary.nsh_max  = ...
                 sp2.boundary(iface).nsh_max + sp3.boundary(iface).nsh_max;
        boundary.nsh      = ...
                 sp2.boundary(iface).nsh + sp3.boundary(iface).nsh;
        boundary.connectivity  = ...
                [sp2.boundary(iface).connectivity; sp3.boundary(iface).connectivity + sp2.boundary(iface).ndof];

        sp2_shape_funs = reshape (sp2.boundary(iface).shape_functions, ...
          [msh.boundary(iface).nqn, sp2.boundary(iface).nsh_max, msh.boundary(iface).nel]);
        sp3_shape_funs = reshape (sp3.boundary(iface).shape_functions, ...
          [msh.boundary(iface).nqn, sp3.boundary(iface).nsh_max, msh.boundary(iface).nel]);

        bnd_shape_functions = zeros (3, msh.boundary(iface).nqn, boundary.nsh_max, msh.boundary(iface).nel);
        bnd_shape_functions(2,:,1:sp2.boundary(iface).nsh_max,:) = sp2_shape_funs;
        bnd_shape_functions(1,:,:,:) = 0;
        bnd_shape_functions(3,:,:,:) = 0;
        bnd_shape_functions(1,:,(sp2.boundary(iface).nsh_max+1):boundary.nsh_max,:) = 0;
        bnd_shape_functions(2,:,(sp2.boundary(iface).nsh_max+1):boundary.nsh_max,:) = 0;
        bnd_shape_functions(3,:,(sp2.boundary(iface).nsh_max+1):boundary.nsh_max,:) = sp3_shape_funs;

        boundary.shape_functions = bnd_shape_functions;

        boundary.comp_dofs{1} = sp1.ndof + sp2.boundary(iface).dofs;
        boundary.comp_dofs{2} = sp1.ndof + sp2.ndof + sp3.boundary(iface).dofs;
        boundary.dofs = [boundary.comp_dofs{:}];

       case {3,4}
        boundary.ncomp    = 3;
        boundary.ndof     = ...
                    sp1.boundary(iface).ndof + sp3.boundary(iface).ndof;
        boundary.ndof_dir = ...
                   [sp1.boundary(iface).ndof_dir; sp3.boundary(iface).ndof_dir];
        boundary.nsh_max  = ...
                    sp1.boundary(iface).nsh_max + sp3.boundary(iface).nsh_max;
        boundary.nsh      = ...
                    sp1.boundary(iface).nsh + sp3.boundary(iface).nsh;
        boundary.connectivity  = ...
                    [sp1.boundary(iface).connectivity; sp3.boundary(iface).connectivity + sp1.boundary(iface).ndof];

        sp1_shape_funs = reshape (sp1.boundary(iface).shape_functions, ...
          [msh.boundary(iface).nqn, sp1.boundary(iface).nsh_max, msh.boundary(iface).nel]);
        sp3_shape_funs = reshape (sp3.boundary(iface).shape_functions, ...
          [msh.boundary(iface).nqn, sp3.boundary(iface).nsh_max, msh.boundary(iface).nel]);

        bnd_shape_functions = zeros (3, msh.boundary(iface).nqn, boundary.nsh_max, msh.boundary(iface).nel);
        bnd_shape_functions(1,:,1:sp1.boundary(iface).nsh_max,:) = sp1_shape_funs;
        bnd_shape_functions(2,:,:,:) = 0;
        bnd_shape_functions(3,:,:,:) = 0;
        bnd_shape_functions(1,:,(sp1.boundary(iface).nsh_max+1):boundary.nsh_max,:) = 0;
        bnd_shape_functions(2,:,(sp1.boundary(iface).nsh_max+1):boundary.nsh_max,:) = 0;
        bnd_shape_functions(3,:,(sp1.boundary(iface).nsh_max+1):boundary.nsh_max,:) = sp3_shape_funs;

        boundary.shape_functions = bnd_shape_functions;

        boundary.comp_dofs{1} = sp1.boundary(iface).dofs;
        boundary.comp_dofs{2} = sp1.ndof + sp2.ndof + sp3.boundary(iface).dofs;
        boundary.dofs = [boundary.comp_dofs{:}];

       case {5,6}
        boundary.ncomp    = 3;
        boundary.ndof     = ...
                    sp1.boundary(iface).ndof + sp2.boundary(iface).ndof;
        boundary.ndof_dir = ...
                   [sp1.boundary(iface).ndof_dir; sp2.boundary(iface).ndof_dir];
        boundary.nsh_max  = ...
                    sp1.boundary(iface).nsh_max + sp2.boundary(iface).nsh_max;
        boundary.nsh      = ...
                    sp1.boundary(iface).nsh + sp2.boundary(iface).nsh;
        boundary.connectivity  = ...
                    [sp1.boundary(iface).connectivity; sp2.boundary(iface).connectivity + sp1.boundary(iface).ndof];

        sp1_shape_funs = reshape (sp1.boundary(iface).shape_functions, ...
          [msh.boundary(iface).nqn, sp1.boundary(iface).nsh_max, msh.boundary(iface).nel]);
        sp2_shape_funs = reshape (sp2.boundary(iface).shape_functions, ...
          [msh.boundary(iface).nqn, sp2.boundary(iface).nsh_max, msh.boundary(iface).nel]);

        bnd_shape_functions = zeros (3, msh.boundary(iface).nqn, boundary.nsh_max, msh.boundary(iface).nel);
        bnd_shape_functions(1,:,1:sp1.boundary(iface).nsh_max,:) = sp1_shape_funs;
        bnd_shape_functions(2,:,:,:) = 0;
        bnd_shape_functions(3,:,:,:) = 0;
        bnd_shape_functions(1,:,(sp1.boundary(iface).nsh_max+1):boundary.nsh_max,:) = 0;
        bnd_shape_functions(2,:,(sp1.boundary(iface).nsh_max+1):boundary.nsh_max,:) = sp2_shape_funs;
        bnd_shape_functions(3,:,(sp1.boundary(iface).nsh_max+1):boundary.nsh_max,:) = 0;

        boundary.shape_functions = bnd_shape_functions;

        boundary.comp_dofs{1} = sp1.boundary(iface).dofs;
        boundary.comp_dofs{2} = sp1.ndof + sp2.boundary(iface).dofs;
        boundary.dofs = [boundary.comp_dofs{:}];
      end

      sp.boundary(iface).ndof_dir = [];
      sp.boundary(iface) = orderfields (boundary, sp.boundary);

      clear bnd_shape_functions
    end
  end

% Map to the physical domain with a curl conserving transform
  sp = sp_curl_transform_3d (sp, msh);

  sp.spfun = @(MSH) sp_bspline_curl_transform_3d ...
                (knots_u1, knots_u2, knots_u3, degree1, degree2, degree3, MSH);

end
