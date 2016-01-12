% SP_VECTOR: Constructor of the class vectorial spaces with a component-wise mapping.
%
%     sp = sp_vector (scalar_spaces, msh, [transform])
%
% INPUTS:
%
%    scalar_spaces: array of space objects, one for each component (see sp_scalar)
%    msh:           mesh object that defines the domain partition and the quadrature rule (see msh_cartesian)
%    transform:     string with the transform to the physical domain, one of 
%                    'grad-preserving' (default), 'curl-preserving' and 'div-preserving'.
%
% OUTPUT:
%
%    sp: object representing the discrete function space of vector-valued functions, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                       DESCRIPTION
%        ncomp           (scalar)                      number of components of the functions of the space (equal to msh.rdim)
%        ncomp_param     (scalar)                      number of components of the functions of the space in the parametric domain (usually equal to msh.ndim)
%        scalar_spaces   (1 x ncomp_param cell array)  array of space objects, one for each component
%        ndof            (scalar)                      total number of degrees of freedom
%        ndof_dir        (ncomp_param x ndim matrix)   for each component, number of degrees of freedom along each direction
%        comp_dofs       (1 x ncomp_param cell array)  indices of the degrees of freedom for each component
%        nsh_max         (scalar)                      maximum number of shape functions per element
%        boundary        (1 x 2*ndim struct array)     struct array representing the space of traces of basis functions on each boundary side
%        transform       (string)                      one of 'grad-preserving', 'curl-preserving' and 'div-preserving'
%        dofs            (1 x ndof vector)             only for boundary spaces, degrees of freedom that do not vanish on the boundary
%        constructor     function handle               function handle to construct the same discrete space in a different msh
%        cumsum_ndof     (1 x ncomp+1 array)           auxiliary array, helpful for indexing
%        cumsum_nsh      (1 x ncomp+1 array)           auxiliary array, helpful for indexing
%
%       METHODS
%       Methods that give a structure with all the functions computed in a certain subset of the mesh
%         sp_evaluate_element_list: compute basis functions (and derivatives) in a given list of elements
%         sp_evaluate_col:          compute basis functions (and derivatives) in one column of the mesh, i.e., fixing the element in the first
%                                    parametric direction
%         sp_precompute:            compute basis functions and derivatives in the whole mesh (memory consuming)       
%         sp_eval_boundary_side:    evaluate the basis functions in one side of the boundary.
%
%       Methods for post-processing, which require a computed vector of degrees of freedom
%         sp_h1_error: compute the error in H1 norm
%         sp_l2_error: compute the error in L2 norm
%         sp_eval:     evaluate the computed solution in a Cartesian grid of points
%         sp_to_vtk:   export the computed solution to a vtk file, using a Cartesian grid of points
%
%       Methods for basic connectivity operations
%         sp_get_basis_functions: compute the functions that do not vanish in a given list of elements
%         sp_get_cells:           compute the cells on which a list of functions do not vanish
%         sp_get_neighbors:       compute the neighbors, functions that share at least one element with a given one
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function sp = sp_vector (scalar_spaces, msh, transform)

  if (nargin == 2)
    transform = 'grad-preserving';
  end

  sp.ncomp = msh.rdim;
  sp.ncomp_param = numel (scalar_spaces);

  switch (transform)
    case {'grad-preserving'}
      if (sp.ncomp ~= sp.ncomp_param)
        error ('sp_vector: the number of components should be equal to the dimension of the physical space')
      end
    case {'curl-preserving'}
      if (numel (scalar_spaces) ~= msh.ndim)
        error ('sp_vector: the dimensions of the space and the mesh do not match')
      end
  end

  sp.scalar_spaces = scalar_spaces;

  sp.nsh_max = sum (cellfun (@(x) x.nsh_max, scalar_spaces));
  sp.ndof    = sum (cellfun (@(x) x.ndof, scalar_spaces));
  sp.ndof_dir = [];
  
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

% Boundary construction
  if (msh.ndim > 1)
    for iside = 1:2*msh.ndim
      for icomp = 1:sp.ncomp_param
        scalar_bnd{icomp} = scalar_spaces{icomp}.boundary(iside);
      end
      
      if (strcmpi (transform, 'grad-preserving'))
        ind = 1:sp.ncomp;
        if (~isempty (msh.boundary))
          sp.boundary(iside) = sp_vector (scalar_bnd(ind), msh.boundary(iside), transform);
        end

      elseif (strcmpi (transform, 'curl-preserving')) % Only tangential components are computed
        ind = setdiff (1:msh.ndim, ceil(iside/2)); % ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
        if (~isempty (msh.boundary))
          sp.boundary(iside) = sp_vector (scalar_bnd(ind), msh.boundary(iside), transform);
        end

      elseif (strcmpi (transform, 'div-preserving')) % Only normal components are computed, and treated as a scalar
        ind = ceil (iside/2); % ind =[1, 1, 2, 2, 3, 3] in 3D, %ind = [1, 1, 2, 2] in 2D;
        if (~isempty (msh.boundary))
          sp_bnd = scalar_bnd{ind};
          if (iside == 1) % This fixes a bug with the use of subsasgn/subsref
            sp.boundary = sp_scalar (sp_bnd.knots, sp_bnd.degree, sp_bnd.weights, msh.boundary(iside), 'integral-preserving');
          else
            sp.boundary(iside) = sp_scalar (sp_bnd.knots, sp_bnd.degree, sp_bnd.weights, msh.boundary(iside), 'integral-preserving');
          end
        end

      else
        error ('sp_vector: unknown transformation')
      end

      dofs = [];
      for icomp = 1:numel(ind)
        new_dofs = sp.cumsum_ndof(ind(icomp)) + scalar_bnd{ind(icomp)}.dofs;
        dofs = union (dofs, new_dofs);
        comp_dofs{icomp} = new_dofs;
      end
      
      sp.boundary(iside).dofs = dofs(:)';
      if (~strcmpi (transform, 'div-preserving'))
        sp.boundary(iside).comp_dofs = comp_dofs;
      end
      
      if (isempty (msh.boundary))
        sp.boundary(iside).ndof = numel (sp.boundary(iside).dofs);
        if (strcmpi (transform, 'curl-preserving')) % Needed for the interfaces on the boundary, in sp_multipatch
          for icomp = 1:numel(ind)
            ndof_dir(icomp,:) = sp.ndof_dir(ind(icomp), ind);
          end
          sp.boundary(iside).ndof_dir = ndof_dir;
        end
      end
    end
    
  elseif (msh.ndim == 1)
    if (strcmpi (transform, 'grad-preserving'))
      sp.boundary(1).dofs = sp.cumsum_ndof(1:end-1)+1;
      sp.boundary(2).dofs = sp.cumsum_ndof(2:end);

      sp.boundary(1).ndof = sp.ncomp_param;
      sp.boundary(2).ndof = sp.ncomp_param;
    else
      sp.boundary = [];
    end
  else
    sp.boundary = [];
  end

  sp.dofs = [];

  sp.transform = transform;

  if (sp.ncomp_param == 2)
    sp.constructor = @(MSH) sp_vector ({scalar_spaces{1}.constructor(MSH), ...
                                        scalar_spaces{2}.constructor(MSH)}, MSH, transform);
  elseif (sp.ncomp_param == 3)
    sp.constructor = @(MSH) sp_vector ({scalar_spaces{1}.constructor(MSH), ...
                                        scalar_spaces{2}.constructor(MSH), ...
                                        scalar_spaces{3}.constructor(MSH)}, MSH, transform);
  elseif (sp.ncomp_param == 1)
    sp.constructor = @(MSH) sp_vector ...
                                      ({scalar_spaces{1}.constructor(MSH)}, MSH, transform);
  end
  
  sp = class (sp, 'sp_vector');


  if (strcmpi (transform, 'div-preserving'))
% Check whether the determinant of the Jacobian is positive
    aux_msh = msh_evaluate_element_list (msh, 1);
    if (any (geopdes_det__ (aux_msh.geo_map_jac) < 0))
      warning ('geopdes:negative_Jacobian', ...
     ['Negative determinant of the Jacobian. This may cause problems with the div-preserving transform. \n' ...
      'For NURBS geometries you can fix this with nrbpermute.'])
    end
  end

end
