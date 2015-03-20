% SP_VECTOR: Constructor of the class vectorial spaces with a component-wise mapping.
%
%     sp = sp_vector (scalar_spaces, msh)
%
% INPUTS:
%
%    scalar_space: array of space objects, one for each component (see sp_bspline or sp_nurbs)
%    msh:          mesh object that defines the domain partition and the quadrature rule (see msh_cartesian)
%
% OUTPUT:
%
%    sp: object representing the discrete function space of vector-valued functions, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                    DESCRIPTION
%        ncomp           (scalar)                  number of components of the functions of the space
%        scalar_spaces   (1 x ncomp cell array)    array of space objects, one for each component
%        ndof            (scalar)                  total number of degrees of freedom
%        ndof_dir        (ncomp x ndim matrix)     for each component, number of degrees of freedom along each direction
%        comp_dofs       (1 x ncomp cell array)    indices of the degrees of freedom for each component
%        nsh_max         (scalar)                  maximum number of shape functions per element
%        boundary        (1 x 2*ndim struct array) struct array representing the space of traces of basis functions on each edge
%
%       METHOD_NAME
%       sp_evaluate_col: compute the basis functions in one column of the mesh (that is, fixing the element in the first parametric direction).
%       sp_eval_boundary_side: evaluate the basis functions in one side of the boundary.
%       sp_precompute:  compute any of the fields related to the discrete
%                       space (except boundary), in all the quadrature points,
%                       as in the space structure from previous versions.
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

function sp = sp_vector (scalar_spaces, msh)

  sp.ncomp = numel (scalar_spaces);
  if (sp.ncomp ~= msh.rdim)
    error ('sp_vector: the number of components should be equal to the dimension of the physical space')
  end
  sp.scalar_spaces = scalar_spaces;

  sp.nsh_max = sum (cellfun (@(x) x.nsh_max, scalar_spaces));
  sp.ndof    = sum (cellfun (@(x) x.ndof, scalar_spaces));
  aux = 0;
  for icomp = 1:sp.ncomp
    sp.comp_dofs{icomp} = aux+(1:scalar_spaces{icomp}.ndof);
    aux = aux + scalar_spaces{icomp}.ndof;
    sp.ndof_dir(icomp,:) = scalar_spaces{icomp}.ndof_dir;
  end
  
  if (~isempty (msh.boundary))
    for iside = 1:numel(msh.boundary)

      for icomp = 1:sp.ncomp
        scalar_bnd{icomp} = scalar_spaces{icomp}.boundary(iside);
      end
      sp.boundary(iside) = sp_vector (scalar_bnd, msh.boundary(iside));
      
      aux = 0;
      dofs = [];
      for icomp = 1:sp.ncomp
%         boundary.ndof_dir(icomp,:) = scalar_bnd{icomp}.ndof_dir;
        new_dofs = aux + scalar_bnd{icomp}.dofs;
        dofs = union (dofs, new_dofs);
        comp_dofs{icomp} = new_dofs;
        aux = aux + scalar_spaces{icomp}.ndof;
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
  sp.shape_function_gradients = [];
  sp.shape_function_divs = [];
  sp.shape_function_curls = [];
  sp.dofs = [];

  if (sp.ncomp == 2)
    sp.constructor = @(MSH) sp_vector ({scalar_spaces{1}.constructor(MSH), ...
                                        scalar_spaces{2}.constructor(MSH)}, MSH);
  elseif (sp.ncomp == 3)
    sp.constructor = @(MSH) sp_vector ({scalar_spaces{1}.constructor(MSH), ...
                                        scalar_spaces{2}.constructor(MSH), ...
                                        scalar_spaces{3}.constructor(MSH)}, MSH);
  end
  sp = class (sp, 'sp_vector');

end
