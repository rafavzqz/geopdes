% SP_VECTOR_CURL_TRANSFORM: Constructor of the class of vectorial spaces with a curl transforming mapping.
%
%     sp = sp_vector_curl_transform (scalar_spaces, msh)
%
% INPUTS:
%
%    scalar_spaces: spline space for each component in the parametric domain (see sp_bspline).
%    msh: mesh object that defines the domain partition and the quadrature rule (see msh_cartesian)
%
% OUTPUT:
%
%    sp: object representing the discrete function space of vector-valued functions, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                      DESCRIPTION
%        scalar_spaces   (1 x ndim space object)     space object for each component
%        ndof            (scalar)                    total number of degrees of freedom
%        ndof_dir        (ndim x ndim matrix)        for each component, number of degrees of freedom along each direction
%        comp_dofs       (1 x ndim cell array)       indices of the degrees of freedom for each component
%        nsh_max         (scalar)                    maximum number of shape functions per element
%        ncomp           (scalar)                    number of vector components in physical space (equal to msh.rdim)
%        ncomp_param     (scalar)                    number of vector components in parametric space (equal to msh.ndim)
%        boundary        (1 x 2*ndim struct array)   struct array representing the space of tangential traces of basis functions on each edge
%
%       METHOD_NAME
%       sp_evaluate_col: compute the basis functions in one column of the mesh (that is, fixing the element in the first parametric direction).
%       sp_eval_boundary_side: evaluate the basis functions in one side of the boundary.
%       sp_precompute:  compute any of the fields related to the discrete
%                       space (except boundary), in all the quadrature points,
%                       as in the space structure from previous versions.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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

function sp = sp_vector_curl_transform (scalar_spaces, msh)

  sp.ncomp = msh.rdim;
  sp.ncomp_param = numel (scalar_spaces);
  if (numel (scalar_spaces) ~= msh.ndim)
    error ('sp_vector_curl_transform: the dimensions of the space and the mesh do not match')
  end
  sp.scalar_spaces = scalar_spaces;

  sp.ndof    = sum (cellfun (@(x) x.ndof, scalar_spaces));
  sp.nsh_max = sum (cellfun (@(x) x.nsh_max, scalar_spaces));
  
  sp.cumsum_ndof(1) = 0;
  sp.cumsum_ndof(2:sp.ncomp_param+1) = cumsum (cellfun (@(x) x.ndof, scalar_spaces));
  sp.cumsum_nsh(1) = 0;
  sp.cumsum_nsh(2:sp.ncomp_param+1) = cumsum (cellfun (@(x) x.nsh_max, scalar_spaces));
  
  aux = 0;
  for icomp = 1:sp.ncomp_param
    sp.comp_dofs{icomp} = aux+(1:scalar_spaces{icomp}.ndof);
    aux = aux + scalar_spaces{icomp}.ndof;
    sp.ndof_dir(icomp,:) = scalar_spaces{icomp}.ndof_dir;
  end
  
  sp.dofs = [];
  if (~isempty (msh.boundary))
    for iside = 1:numel(msh.boundary)
      for icomp = 1:sp.ncomp_param
        scalar_bnd{icomp} = scalar_spaces{icomp}.boundary(iside);
      end

%%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
      ind = setdiff (1:msh.ndim, ceil(iside/2)); 
      sp.boundary(iside) = sp_vector_curl_transform (scalar_bnd(ind), msh.boundary(iside));
      
      dofs = [];
      for icomp = 1:numel(ind)
        new_dofs = sp.cumsum_ndof(ind(icomp)) + scalar_bnd{ind(icomp)}.dofs;
        dofs = union (dofs, new_dofs);
        comp_dofs{icomp} = new_dofs;
      end
      sp.boundary(iside).dofs = dofs(:)';
      sp.boundary(iside).comp_dofs = comp_dofs;
    end
  else
    sp.boundary = [];
  end

  sp.nsh = [];
  sp.connectivity = [];
  sp.shape_functions = [];
  sp.shape_function_curls = [];

  if (sp.ncomp_param == 2)
    sp.constructor = @(MSH) sp_vector_curl_transform ...
                                      ({scalar_spaces{1}.constructor(MSH), ...
                                        scalar_spaces{2}.constructor(MSH)}, MSH);
  elseif (sp.ncomp_param == 3)
    sp.constructor = @(MSH) sp_vector_curl_transform ...
                                      ({scalar_spaces{1}.constructor(MSH), ...
                                        scalar_spaces{2}.constructor(MSH), ...
                                        scalar_spaces{3}.constructor(MSH)}, MSH);
  elseif (sp.ncomp_param == 1)
    sp.constructor = @(MSH) sp_vector_curl_transform ...
                                      ({scalar_spaces{1}.constructor(MSH)}, MSH);
  end
  sp = class (sp, 'sp_vector_curl_transform');

end
