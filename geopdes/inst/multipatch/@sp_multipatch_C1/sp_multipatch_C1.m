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
%        ncomp           (scalar)                      number of components of the functions of the space (always equal to one)
%        ndof            (scalar)                      total number of degrees of freedom after gluing patches together
%        ndof_per_patch  (1 x npatch array)            number of degrees of freedom per patch, without gluing
% interior_dofs_per_patch
% ndof_interior
% ndof_interface
%        sp_patch        (1 x npatch cell-array)       the input spaces, one space object for each patch (see sp_scalar and sp_vector)
%        gnum            (1 x npatch cell-array)       global numbering of the degress of freedom (see mp_interface)
%        constructor     function handle               function handle to construct the same discrete space in a different msh
%
%       METHODS
%       Methods that give a structure with all the functions computed in a certain subset of the mesh
% % % %         sp_evaluate_element_list: compute basis functions (and derivatives) in a given list of elements
%
%       Methods for post-processing, that require a computed vector of degrees of freedom
%         sp_l2_error:    compute the error in L2 norm
%         sp_h1_error:    compute the error in H1 norm
%         sp_h2_error:    compute the error in H2 norm
%         sp_to_vtk:      export the computed solution to a pvd file, using a Cartesian grid of points on each patch
%
%       Methods for basic connectivity operations
%         sp_get_basis_functions: compute the functions that do not vanish in a given list of elements
%         sp_get_cells:           compute the cells on which a list of functions do not vanish
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

% XXX SHOULD ADD A CHECK ABOUT DEGREE AND REGULARITY

  if (~all (cellfun (@(x) isa (x, 'sp_scalar'), spaces)))
    error ('All the spaces in the array should be of the same class')
  end
  aux = struct ([spaces{:}]);

  sp.npatch = numel (spaces);
  if (sp.npatch ~= msh.npatch)
    error ('The list of spaces does not correspond to the mesh')
  end

% XXX FIX THIS, TO COMPUTE ALSO THE BOUNDARY
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
% XXX FOR NOW ONLY FOR TWO PATCHES, AND WITH THE ORDERING AS IN THE PAPER
% XXX THE TWO PATCHES ARE ALSO ORIENTED AS IN THE PAPER
  sp.ndof_interior = 0;
  for iptc = 1:sp.npatch
    ndof_dir = sp.sp_patch{iptc}.ndof_dir;
    xx = 3:ndof_dir(1); yy = 1:ndof_dir(2);
    [XX, YY] = ndgrid (xx, yy);
    interior_dofs = sub2ind (ndof_dir, XX, YY);
    sp.interior_dofs_per_patch{iptc} = interior_dofs; % An array with the indices
    sp.ndof_interior = sp.ndof_interior + numel (interior_dofs);
  end
  
  [ndof, CC] = compute_coefficients (sp, msh, geometry, interfaces);
  
  sp.ndof_interface = ndof; % The number of functions in V^2
  sp.ndof = sp.ndof_interior + sp.ndof_interface;


% Computation of the coefficients for basis change
% The matrix Cpatch{iptc} is a matrix of size ndof_per_patch(iptc) x ndof
% The coefficients for basis change have been stored in CC

  Cpatch = cell (sp.npatch, 1);
  numel_interior_dofs = cellfun (@numel, sp.interior_dofs_per_patch);
  for iptc = 1:sp.npatch
    Cpatch{iptc} = sparse (sp.ndof_per_patch(iptc), sp.ndof);
    global_indices = sum (numel_interior_dofs(1:iptc-1)) + (1:numel_interior_dofs(iptc));
    Cpatch{iptc}(sp.interior_dofs_per_patch{iptc}, global_indices) = ...
      speye (numel (sp.interior_dofs_per_patch{iptc}));

    global_indices = (sp.ndof_interior+1):sp.ndof;
    Cpatch{iptc}(:,global_indices) = CC{iptc};
  end  

  sp.interfaces = interfaces;
  sp.Cpatch = Cpatch;

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


function [ndof, CC] = compute_coefficients (space, msh, geometry, interfaces)

  for iref = 1:numel(interfaces)
    patch(1) = interfaces(iref).patch1;
    patch(2) = interfaces(iref).patch2;
    side(1) = interfaces(iref).side1;
    side(2) = interfaces(iref).side2;

    
% The knot vectors for the N0 and N1 basis functions, as in Mario's notation
% I assume that the regularity is degree-2, otherwise things become more complicated
% Only univariate knot vectors are computed
    if (side(1) < 3)
      knots = space.sp_patch{patch(1)}.knots{2};
      degree = space.sp_patch{patch(1)}.degree(2);
      nel_univ = msh.msh_patch{patch(1)}.nel_dir(2);
      degu = space.sp_patch{patch(1)}.degree(1);
      knt = unique (space.sp_patch{patch(1)}.knots{1});
      tau1 = knt(2) - knt(1);
    else
      knots = space.sp_patch{patch(1)}.knots{1};
      degree = space.sp_patch{patch(1)}.degree(1);
      nel_univ = msh.msh_patch{patch(1)}.nel_dir(1);
      degu = space.sp_patch{patch(1)}.degree(2);
      knt = unique (space.sp_patch{patch(1)}.knots{2});
      tau1 = knt(2) - knt(1);
    end
    regularity = degree - 2;
    
    nel_geo = numel (unique (geometry(patch(1)).boundary(side(1)).nurbs.knots)) - 1;
    nsub = nel_univ / nel_geo;
    
    knots0 = kntrefine (geometry(patch(1)).boundary(side(1)).nurbs.knots, nsub-1, degree, regularity+1);
    knots1 = kntrefine (geometry(patch(1)).boundary(side(1)).nurbs.knots, nsub-1, degree-1, regularity);

% Compute the Greville points, and the auxiliary mesh and space objects
    for ii = 1:2 % The two patches (L-R)
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

% For now I assume that the orientation is as in the paper, and we do not need any reordering
%XXXX    [sp_bnd(2), msh_side(2)] = reorder_elements_and_quad_points (sp_bnd(2), msh_side(2), interfaces(iref), ndim);
% msh_side contains only the univariate parametrization of the boundary (dependence on u)
% msh_side_int contains information for the bivariate parametrization (dependence on u and v)
% sp_grev contains the value and derivatives of the basis functions, at the Greville points
      msh_side(ii) = msh_eval_boundary_side (msh_grev, side(ii));
      msh_side_int{ii} = msh_boundary_side_from_interior (msh_grev, side(ii));
      sp_aux = space.sp_patch{patch(ii)}.constructor (msh_side_int{ii});
      msh_side_int{ii} = msh_precompute (msh_side_int{ii});
      sp_grev(ii) = struct (sp_precompute_param (sp_aux, msh_side_int{ii}, 'value', true, 'gradient', true));

% The univariate spaces for the N0 and N1 basis functions      
      sp0 = sp_bspline (knots0, degree, msh_grev.boundary(side(ii)));
      sp1 = sp_bspline (knots1, degree-1, msh_grev.boundary(side(ii)));
      sp0_struct = sp_precompute_param (sp0, msh_grev.boundary(side(ii)), 'value', true, 'gradient', true);
      sp1_struct = sp_precompute_param (sp1, msh_grev.boundary(side(ii)), 'value', true, 'gradient', true);
      ind = [2 2 1 1];
      knotsn = sp_aux.knots{ind(side(ii))};
      spn = sp_bspline (knotsn, degree, msh_grev.boundary(side(ii)));
      spn_struct = sp_precompute_param (spn, msh_grev.boundary(side(ii)), 'value', true, 'gradient', true);
% Matrix for the linear systems, (14)-(16) in Mario's notes
      A = sparse (msh_side(ii).nel, msh_side(ii).nel);
      for jj = 1:msh_side(ii).nel
        A(jj,spn_struct.connectivity(:,jj)) = squeeze (spn_struct.shape_functions(:,:,jj));
      end

      alpha{ii} = geopdes_det__ (msh_side_int{ii}.geo_map_jac);
      numerator = reshape (sum (msh_side_int{ii}.geo_map_jac(:,1,:,:) .* msh_side_int{ii}.geo_map_jac(:,2,:,:), 1), msh_side(ii).nel, 1);
      denominator = reshape (sum (msh_side_int{ii}.geo_map_jac(:,2,:,:) .* msh_side_int{ii}.geo_map_jac(:,2,:,:), 1), msh_side(ii).nel, 1);
      beta{ii} = numerator ./ denominator;

% RHS for the first linear system, (14) in Mario's notes
      rhss = sparse (msh_side(ii).nel, sp0_struct.ndof);
      for jj = 1:msh_side(ii).nel
        rhss(jj,sp0_struct.connectivity(:,jj)) = squeeze (sp0_struct.shape_functions(:,:,jj));
      end
      coeff0{ii} = A \ rhss;
      coeff0{ii}(abs(coeff0{ii}) < 1e-12) = 0; % Make more sparse

% RHS for the second linear system, (15) in Mario's notes
      rhsb = sparse (msh_side(ii).nel, sp0_struct.ndof);
      val = squeeze (spn_struct.shape_function_gradients(:,:,2,1)); % XXX THIS IS ONLY VALID FOR THE ORDERING USED IN THE PAPER
      val = val * (tau1 / degu)^2;
      for jj = 1:msh_side(ii).nel
        val_aux = val * beta{ii}(jj);
        rhsb(jj,sp0_struct.connectivity(:,jj)) = squeeze (sp0_struct.shape_function_gradients(:,:,:,jj)) * val_aux;
      end
      rhsb = rhsb + rhss;
      coeff1{ii} = A \ rhsb;
      coeff1{ii}(abs(coeff1{ii}) < 1e-12) = 0; % Make more sparse
      
% RHS for the third linear system, (16) in Mario's notes
      rhsc = sparse (msh_side(ii).nel, sp1_struct.ndof);
      val = squeeze (spn_struct.shape_function_gradients(:,:,2,1)); % XXX THIS IS ONLY VALID FOR THE ORDERING USED IN THE PAPER
      val = val * (tau1 / degu)^2;
      for jj = 1:msh_side(ii).nel
        val_aux = val * alpha{ii}(jj);
        rhsc(jj,sp1_struct.connectivity(:,jj)) = squeeze (sp1_struct.shape_functions(:,:,jj)) * val_aux;
      end
      coeff2{ii} = A \ rhsc;
      coeff2{ii}(abs(coeff2{ii}) < 1e-12) = 0; % Make more sparse
      
% Pass the coefficients to the tensor product basis
% The numbering (ndof) only works for the two patch case, for now
      ndof = sp0_struct.ndof + sp1_struct.ndof;
      CC{ii} = sparse (space.ndof_per_patch(patch(ii)), ndof);
      
      ind0 = sub2ind (space.sp_patch{patch(ii)}.ndof_dir, ones(1,spn.ndof), 1:spn.ndof);
      ind1 = sub2ind (space.sp_patch{patch(ii)}.ndof_dir, 2*ones(1,spn.ndof), 1:spn.ndof);
%       ind2 = sub2ind (space.sp_patch{patch(ii)}.ndof_dir, 2*ones(1,spn.ndof), 1:spn.ndof);

      CC{ii}(ind0,1:sp0_struct.ndof) = coeff0{ii};
      CC{ii}(ind1,1:sp0_struct.ndof) = coeff1{ii};
      CC{ii}(ind1,sp0_struct.ndof+1:ndof) = coeff2{ii};
    end


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
