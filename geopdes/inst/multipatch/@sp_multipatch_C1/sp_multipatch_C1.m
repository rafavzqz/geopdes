% SP_MULTIPATCH_C1: Constructor of the class for multipatch spaces with C1 continuity
%
% BETA VERSION. For now, it will only work with two patches
%
%     sp = sp_multipatch_C1 (spaces, msh, interfaces)
%%%XXXX     sp = sp_multipatch_C1 (spaces, msh, interfaces, boundary_interfaces)
%
% INPUTS:
%
%    spaces:     cell-array of space objects, one for each patch (see sp_scalar, sp_vector)
%    msh:        mesh object that defines the multipatch mesh (see msh_multipatch)
%    interfaces: information of connectivity between patches (see mp_geo_load)
%%%XXXX    boundary_interfaces: information of connectivity between boundary patches (see mp_geo_load)
%
% OUTPUT:
%
%    sp: object representing the discrete function space of vector-valued functions, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                       DESCRIPTION
%        npatch          (scalar)                      number of patches
%XXXX ncomp = 1        ncomp           (scalar)                      number of components of the functions of the space (equal to msh.rdim)
%        ndof            (scalar)                      total number of degrees of freedom after gluing patches together
%        ndof_per_patch  (1 x npatch array)            number of degrees of freedom per patch, without gluing
% interior_dofs_per_patch
% ndof_interior
% ndof_interface
%        sp_patch        (1 x npatch cell-array)       the input spaces, one space object for each patch (see sp_scalar and sp_vector)
%XXXX        transform       (string)                      one of 'grad-preserving', 'curl-preserving' and 'div-preserving'
%        gnum            (1 x npatch cell-array)       global numbering of the degress of freedom (see mp_interface)
%XXXX        dofs_ornt       (1 x npatch cell-array)       global orientation of the degrees for freedom, for curl-conforming and div-conforming spaces
%XXXX        boundary        (1 x 1 object)                a (ndim-1) dimensional "sp_multipatch" object for the whole boundary 
%XXXX        dofs            (1 x ndof vector)             only for boundary spaces, degrees of freedom that do not vanish on the boundary
%XXXX        orientation     (1 x ndof vector)             only for boundary spaces with 'curl-preserving' transform, global orientation
%XXXX                                                       of the boundary dofs with respect to the volumetric one.
%        constructor     function handle               function handle to construct the same discrete space in a different msh
%
% % % %       METHODS
% % % %       Methods that give a structure with all the functions computed in a certain subset of the mesh
% % % %         sp_evaluate_element_list: compute basis functions (and derivatives) in a given list of elements
% % % %
% % % %       Methods for post-processing, that require a computed vector of degrees of freedom
% % % %         sp_h1_error:    compute the error in H1 norm
% % % %         sp_l2_error:    compute the error in L2 norm
% % % %         sp_hcurl_error: compute the error in H(curl) norm
% % % %         sp_to_vtk:      export the computed solution to a pvd file, using a Cartesian grid of points on each patch
% % % %
% % % %       Methods for basic connectivity operations
% % % %         sp_get_basis_functions: compute the functions that do not vanish in a given list of elements
% % % %         sp_get_cells:           compute the cells on which a list of functions do not vanish
% % % %         sp_get_neighbors:       compute the neighbors, functions that share at least one element with a given one
%
% Copyright (C) 2015, 2017 Rafael Vazquez
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

function sp = sp_multipatch_C1 (spaces, msh, geometry, interfaces, boundary_interfaces)

  if (~all (cellfun (@(x) isa (x, 'sp_scalar'), spaces)))
    error ('All the spaces in the array should be of the same class')
  end
  aux = struct ([spaces{:}]);

  sp.npatch = numel (spaces);
  if (sp.npatch ~= msh.npatch)
    error ('The list of spaces does not correspond to the mesh')
  end

% XXXX FIX THIS, TO COMPUTE ALSO THE BOUNDARY
  if (msh.ndim ~= 2 || msh.rdim ~= 2)
    error ('Only implemented for planar surfaces')
  end


  sp.ncomp = spaces{1}.ncomp;
  sp.transform = spaces{1}.transform;
  
%   if (~all ([aux.ncomp] == sp.ncomp))
  if (~all ([aux.ncomp] == 1))
    error ('The number of components should be the same for all the spaces, and equal to one')  
  end
  for iptc = 1:sp.npatch
%     if (~strcmpi (spaces{iptc}.transform, sp.transform))
    if (~strcmpi (spaces{iptc}.transform, 'grad-preserving'))
      error ('The transform to the physical domain should be the same for all the spaces, and the grad-preserving one')
    end
    if (~strcmpi (spaces{iptc}.space_type, 'spline'))
      error ('C1 continuity is only implemented for splines, not for NURBS')
    end
  end
  
  sp.ndof = 0;
  sp.ndof_per_patch = [aux.ndof];
  sp.sp_patch = spaces;

% Computation of the number of degrees of freedom
% I need to give a global numbering to the C^1 basis functions
% I start numbering those away from the interface (V^1) patch by patch
% And then generate the numbering for the functions close to the interface (V^2)

% Compute the local indices of the functions in V^1
%  and sum them up to get the whole space V^1
  sp.ndof_interior = 0;
  for iptc = 1:sp.npatch
    sp.interior_dofs_per_patch{iptc} = []; % An array with the indices
    sp.ndof_interior = sp.ndof_interior + numel (sp.interior_dofs_per_patch{iptc});
  end
  sp.ndof_interface = 0; % The number of functions in V^2
  sp.ndof = sp.ndof_interior + sp.ndof_interface;


% Computation of the coefficients for basis change
% The matrix Cpatch{iptc} is a matrix of size ndof_per_patch(iptc) x ndof
% Here we should call the function compute_coefficients below


  Cpatch = cell (sp.npatch, 1);
  numel_interior_dofs = cellfun (@numel, sp.interior_dofs_per_patch);
  for iptc = 1:sp.npatch
    Cpatch{iptc} = sparse (sp.ndof_per_patch(iptc), sp.ndof);
    global_indices = sum (numel_interior_dofs(1:iptc-1)) + (1:numel_interior_dofs(iptc));
    Cpatch{iptc}(sp.interior_dofs_per_patch{iptc}, global_indices) = ...
      speye (numel (sp.interior_dofs_per_patch{iptc}));

    global_indices = (sp.ndof_interior+1):sp.ndof;
    Cpatch{iptc}(:,global_indices) = []; % The coefficients computed by Mario

disp ('We should fill Cpatch using the coefficients, computed as in Section 7')
disp('Mario, use F11 (step In) or F10 (Step) to see where we are')
keyboard

    compute_coefficients (sp, msh, geometry, interfaces)

  end  

  sp.interfaces = interfaces;

% XXXX No boundary for now
% % Boundary construction
%   if (nargin == 4 && ~isempty (msh.boundary) && ~isempty (spaces{1}.boundary))
%     sp_bnd = cell (msh.boundary.npatch, 1);
%     for iptc = 1:msh.boundary.npatch
%       patch_number = msh.boundary.patch_numbers(iptc);
%       side_number  = msh.boundary.side_numbers(iptc);
%       sp_bnd{iptc} = spaces{patch_number}.boundary(side_number);
%     end
%     sp.boundary = sp_multipatch (sp_bnd, msh.boundary, boundary_interfaces);
%     
%     dofs = zeros (sp.boundary.ndof, 1);
%     boundary_orientation = []; %boundary_orient = zeros (sp.boundary.ndof, 1);
%     for iptc = 1:msh.boundary.npatch
%       patch_number = msh.boundary.patch_numbers(iptc);
%       side_number  = msh.boundary.side_numbers(iptc);
%       dofs(sp.boundary.gnum{iptc}) = sp.gnum{patch_number}(sp.sp_patch{patch_number}.boundary(side_number).dofs);
%       if (~isempty (sp.boundary.dofs_ornt))
%         boundary_orientation(sp.boundary.gnum{iptc}) = sp.boundary.dofs_ornt{iptc} .* ...
%           sp.dofs_ornt{patch_number}(sp.sp_patch{patch_number}.boundary(side_number).dofs);
%       end
%     end
%     sp.boundary.dofs = dofs;
%     sp.boundary.boundary_orientation = boundary_orientation;
%     
%   else
%     sp.boundary = [];
%   end
%   
%   sp.dofs = [];
%   sp.boundary_orientation = [];
  
  sp.constructor = @(MSH) sp_multipatch_C1 (patches_constructor(spaces, MSH), MSH, geometry, interfaces);
    function spaux = patches_constructor (spaces, MSH)
      for ipatch = 1:MSH.npatch
        spaux{ipatch} = spaces{ipatch}.constructor(MSH.msh_patch{ipatch});
      end
    end

  sp = class (sp, 'sp_multipatch_C1');
  
end


function compute_coefficients (space, msh, geometry, interfaces)

  for iref = 1:numel(interfaces)
    patch(1) = interfaces(iref).patch1;
    patch(2) = interfaces(iref).patch2;
    side(1) = interfaces(iref).side1;
    side(2) = interfaces(iref).side2;

% Compute the Greville points, and the auxiliary mesh and space objects
    for ii = 1:2 % The two patches
      brk = cell (1,msh.ndim);
      knots = space.sp_patch{patch(ii)}.knots;
      for idim = 1:msh.ndim
        grev_pts{idim} = aveknt (knots{idim}, space.sp_patch{patch(ii)}.degree(idim)+1); 
        grev_pts{idim} = grev_pts{idim}(:)';
        if (numel(grev_pts{idim}) > 1)
          brk{idim} = [knots{idim}(1), grev_pts{idim}(1:end-1) + diff(grev_pts{idim})/2, knots{idim}(end)];
        else
          brk{idim} = [knots{idim}(1) knots{idim}(end)];
        end
      end
      msh_grev = msh_cartesian (brk, grev_pts, [], geometry(patch(ii)), 'boundary', true, 'der2',false);

% msh_side contains only the univariate parametrization of the boundary (dependence on u)
% msh_side_int contains information for the bivariate parametrization (dependence on u and v)
% sp_grev contains the value and derivatives of the basis functions, at the Greville points
      msh_side(ii) = msh_eval_boundary_side (msh_grev, side(ii));
      msh_side_int{ii} = msh_boundary_side_from_interior (msh_grev, side(ii));
      sp_aux = space.sp_patch{patch(ii)}.constructor (msh_side_int{ii});
      msh_side_int{ii} = msh_precompute (msh_side_int{ii});
      sp_grev(ii) = struct (sp_precompute_param (sp_aux, msh_side_int{ii}, 'value', true, 'gradient', true));
    end
% For now I assume that the orientation is as in the paper, and we do not need any reordering
%XXXX    [sp_bnd(2), msh_side(2)] = reorder_elements_and_quad_points (sp_bnd(2), msh_side(2), interfaces(iref), ndim);

  end
  
  
%%%%%%%%%% FOR MARIO, TRYING TO SUMMARIZE %%%%%%%%%%%%%%%%
% In the inner loop (for ii = 1:2), we consider the two patches. ii = L,R
% msh_side_int{ii} is a structure that contains the information about the parametrization F^(ii)
%  In particular, it contains the value (geo_map), the derivatives (geo_map_jac) 
%  and the jacobian (jacdet) at the Greville points. The last index
%  indicates the Greville point.
% sp_grev(ii) is a structure with the value of the basis functions
% (.shape_functions) at the Greville points. The last index indicates the Greville point.
% The previous one indicates the basis functions that are non-null on the
% element. This also includes functions from the interior.
% The number of the non-vanishing basis functions is in connectivity

end
