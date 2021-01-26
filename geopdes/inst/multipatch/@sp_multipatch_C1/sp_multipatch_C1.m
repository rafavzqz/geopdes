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
%    vertices: information about interfaces and patches neighbouring each vertex (to be implemented)
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
%         sp_evaluate_element_list: compute basis functions (and derivatives) in a given list of elements
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
%         sp_get_neighbors:       compute the neighbors, functions that share at least one element with a given one
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

function sp = sp_multipatch_C1 (spaces, msh, geometry, interfaces, boundaries)

  if (~all (cellfun (@(x) isa (x, 'sp_scalar'), spaces)))
    error ('All the spaces in the array should be of the same class')
  end
  aux = struct ([spaces{:}]);

  sp.npatch = numel (spaces);
  if (sp.npatch ~= msh.npatch)
    error ('The list of spaces does not correspond to the mesh')
  end

  if (msh.ndim ~= 2 || msh.rdim ~= 2)
    error ('Only implemented for planar surfaces')
  end
  

   for iptc = 1:numel(geometry)
%     if (any (geometry(iptc).nurbs.order > 2))
%       error ('For now, only bilinear patches are implemented')
%     end
    knots = spaces{iptc}.knots;
    breaks = cellfun (@unique, knots, 'UniformOutput', false);
    mult = cellfun (@histc, knots, breaks, 'UniformOutput', false);
    if (any ([mult{:}] < 2))
      error ('The regularity should be at most degree minus two')
    end
    for idim = 1:2
      if (any (mult{idim}(2:end-1) > spaces{iptc}.degree(idim) - 1))
        error ('The regularity should not be lower than one')
      end
    end
  end
  

  sp.ncomp = spaces{1}.ncomp;
  sp.transform = spaces{1}.transform;
  
  if (~all ([aux.ncomp] == 1))
    error ('The number of components should be the same for all the spaces, and equal to one')  
  end
  for iptc = 1:sp.npatch
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
  
% Assuming that the starting space has degree p and regularity r, 
% r <= p-2, we compute the knot vectors of the auxiliary spaces:
%  knots0: degree p, regularity r+1.
%  knots1: degree p-1, regularity r.

  knots0 = cell (sp.npatch, 1); knots1 = knots0;
  for iptc = 1:sp.npatch
    knots = spaces{iptc}.knots;
    breaks = cellfun (@unique, knots, 'UniformOutput', false);
    for idim = 1:msh.ndim
      mult = histc (knots{idim}, breaks{idim});
      mult0{idim} = mult; mult0{idim}(2:end-1) = mult(2:end-1) - 1;
      mult1{idim} = mult - 1;
    end
    knots0{iptc} = kntbrkdegmult (breaks, spaces{iptc}.degree, mult0);
    knots1{iptc} = kntbrkdegmult (breaks, spaces{iptc}.degree-1, mult1);
  end
  sp.knots0_patches = knots0;
  sp.knots1_patches = knots1;

% Computation of the number of degrees of freedom
% We need to give a global numbering to the C^1 basis functions
% We start numbering those away from the interface (V^1) patch by patch
% And then generate the numbering for the functions close to the interface (V^2)

% Compute the local indices of the functions in V^1
%  and sum them up to get the whole space V^1

% if (numel (interfaces) > 1 || msh.npatch > 2)
%   error ('For now, the implementation only works for two patches')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I HAVE CHANGED STARTING FROM HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We could store the info of interfaces per patch, as in mp_interface

  sp.ndof_interior = 0;
  [interfaces_all, ~] = vertices_struct(boundaries, interfaces);
  for iptc = 1:sp.npatch
    interior_dofs = 1:spaces{iptc}.ndof;
    for intrfc = 1:numel(interfaces_all)
      patches = [interfaces_all(intrfc).patch1, interfaces_all(intrfc).patch2];
      sides = [interfaces_all(intrfc).side1, interfaces_all(intrfc).side2];
      [is_interface,position] = ismember (iptc, patches);
      if (is_interface)
        sp_bnd = spaces{iptc}.boundary(sides(position));
        interior_dofs = setdiff (interior_dofs, [sp_bnd.dofs, sp_bnd.adjacent_dofs]);
      end
    end
    sp.interior_dofs_per_patch{iptc} = interior_dofs; % An array with the indices
    sp.ndof_interior = sp.ndof_interior + numel (interior_dofs);
  end
  
% EXPECTED OUTPUT FROM compute_coefficients  
% ndof_per_interface: number of edge functions on each interface, array of
%    size 1 x numel(interfaces);
% CC_edges: cell array of size 2 x numel(interfaces), the two corresponds
%   to the two patches on the interface. The matrix CC_edges{ii,jj} has size
%      sp.ndof_per_patch(patch) x ndof_per_interface(jj)
%      with patch = interfaces(jj).patches(ii);
% ndof_per_vertex: number of vertex functions on each vertex. An array of size numel(vertices)
% CC_vertices: cell array of size npatch x numel(vertices)
%    The matrix CC_vertices{ii,jj} has size
%    sp.ndof_per_patch(patch) x ndof_per_vertex{jj}
%      with patch being the index of the ii-th patch containing vertex jj;
%
  
%   [ndof, CC] = compute_coefficients (sp, msh, geometry, interfaces);  
  [ndof_per_interface, CC_edges, ndof_per_vertex, CC_vertices] = ...
    compute_coefficients (sp, msh, geometry, interfaces, boundaries);
%keyboard
  
  sp.ndof_edges = sum(ndof_per_interface); % Total number of edge functions
  sp.ndof_vertices = sum (ndof_per_vertex); % Total number of vertex functions
  sp.ndof = sp.ndof_interior + sp.ndof_edges + sp.ndof_vertices;

% Computation of the coefficients for basis change
% The matrix Cpatch{iptc} is a matrix of size ndof_per_patch(iptc) x ndof
% The coefficients for basis change have been stored in CC_*

  Cpatch = cell (sp.npatch, 1);
  numel_interior_dofs = cellfun (@numel, sp.interior_dofs_per_patch);
  for iptc = 1:sp.npatch
    Cpatch{iptc} = sparse (sp.ndof_per_patch(iptc), sp.ndof);
    global_indices = sum (numel_interior_dofs(1:iptc-1)) + (1:numel_interior_dofs(iptc));
    Cpatch{iptc}(sp.interior_dofs_per_patch{iptc}, global_indices) = ...
      speye (numel (sp.interior_dofs_per_patch{iptc}));    
  end
  
  for intrfc = 1:numel(interfaces)
    global_indices = sp.ndof_interior + sum(ndof_per_interface(1:intrfc-1)) + (1:ndof_per_interface(intrfc));
    %for iptc_on_interface = 1:2
      iptc1 = interfaces(intrfc).patch1;
      iptc2 = interfaces(intrfc).patch2;
      Cpatch{iptc1}(:,global_indices) = CC_edges{1,intrfc};
      Cpatch{iptc2}(:,global_indices) = CC_edges{2,intrfc};
    %end
  end

% Vertices and patches_on_vertex are not defined yet. For now, this only works for one extraordinary point
% The information of which patches share the vertex can be computed with the help of mp_interface
  for ivrt = 1%:numel(vertices)
    global_indices = sp.ndof_interior + sp.ndof_edges + sum(ndof_per_vertex(1:ivrt-1)) + (1:ndof_per_vertex(ivrt));
    for iptc = 1:sp.npatch %patches_on_vertex (TO BE CHANGED)
      Cpatch{iptc}(:,global_indices) = CC_vertices{iptc,ivrt};
    end
  end

%  sp.patches_on_vertex = patches_on_vertex; TO BE DONE
  sp.interfaces = interfaces;
  sp.Cpatch = Cpatch;
  sp.geometry = geometry; % I store this for simplicity
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I HAVE FINISHED MY CHANGES HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


% There are some issues related to the orientation.
% I am taking the absolute value of alpha
% And for val_grad, I have to change the sign for coeff2, but not for coeff1
% But everything seems to work!!!

function [ndof_per_interface, CC_edges, ndof_per_vertex, CC_vertices] = compute_coefficients (space, msh, geometry, interfaces, boundaries) %based on first Mario's notes (not refinement mask)

[interfaces_all, vertices]=vertices_struct(boundaries, interfaces); 

p = space.sp_patch{1}.degree(1);
k = numel(msh.msh_patch{1}.breaks{1})-2;
n= space.sp_patch{1}.sp_univ(1).ndof;

ndof_per_interface=(2*n-2*k-11)*ones(1,numel(interfaces_all)); %works only if the space is degree and regularity are the same everywhere
ndof_per_vertex=6*ones(1,numel(vertices));

all_alpha0=zeros(numel(interfaces_all),2);
all_alpha1=zeros(numel(interfaces_all),2);
all_beta0=zeros(numel(interfaces_all),2);
all_beta1=zeros(numel(interfaces_all),2);
all_t0=zeros(numel(interfaces_all),2);

%Initialize the cell array of CC_vertices matrices
for ii=1:space.npatch
    for jj=numel(vertices)
        CC_vertices{ii,jj} = sparse (space.ndof_per_patch(ii), ndof_per_vertex(jj));
    end
end
%Initialize the cell array of CC_edges matrices
for jj=1:numel(interfaces_all)
    if isempty(interfaces_all(jj).patch1)
        CC_edges{1,jj} = [];
    else
        CC_edges{1,jj} = sparse (space.ndof_per_patch(interfaces_all(jj).patch1), ndof_per_interface(jj));
    end
    if isempty(interfaces_all(jj).patch2)
        CC_edges{2,jj} = [];
    else
        CC_edges{2,jj} = sparse (space.ndof_per_patch(interfaces_all(jj).patch2), ndof_per_interface(jj));
    end
end

%computing for each patch all the derivatives we possibly need to compute t,d, and sigma
brk = cell (1,msh.ndim);
for j=1:space.npatch
    knots = space.sp_patch{j}.knots;
    for idim = 1:msh.ndim
        brk{idim}=[knots{idim}(1) knots{idim}(end)]; %is this correct?
    end
    %the following points correspond to the four vertices of the patch
    pts{1}=[0 1]';
    pts{2}=[0 1]';%pts{2}=[0 1/2 1]'
    msh_pts_der1 = msh_cartesian (brk, pts, [], geometry(j),'boundary', true, 'der2',false);
    msh_der1{j}=msh_precompute (msh_pts_der1); 
    derivatives1{j}=msh_der1{j}.geo_map_jac; %rdim x ndim x (n_pts{1}x n_pts{2}) (rdim->physical space, ndim->parametric space)
    msh_pts_der2 = msh_cartesian (brk, pts, [], geometry(j),'boundary', true, 'der2',true);
    msh_der2{j}=msh_precompute (msh_pts_der2);     
    derivatives2{j}=msh_der1{j}.geo_map_der2; %rdim x ndim x ndim x n_pts{1}x n_pts{2}
    %In this case (rdim=2)
    %D_u F^j(0,0)=squeeze(derivatives1{j}(:,1,1,1)) %this contains both components in the physical space
    %D_u F^j(0,1)=squeeze(derivatives1{j}(:,1,1,2))
    %D_u F^j(1,0)=squeeze(derivatives1{j}(:,1,2,1))
    %D_u F^j(1,1)=squeeze(derivatives1{j}(:,1,2,2))
    %D_v F^j(0,0)=squeeze(derivatives1{j}(:,2,1,1))
    %D_v F^j(0,1)=squeeze(derivatives1{j}(:,2,1,2))
    %D_v F^j(1,0)=squeeze(derivatives1{j}(:,2,2,1))
    %D_v F^j(1,1)=squeeze(derivatives1{j}(:,2,2,2))    
    
    %D_uu F^j(0,0)=squeeze(derivatives2{j}(:,1,1,1,1)) %this contains both components in the physical space
    %D_uu F^j(0,1)=squeeze(derivatives2{j}(:,1,1,1,2))
    %D_uu F^j(1,0)=squeeze(derivatives2{j}(:,1,1,2,1))
    %D_uu F^j(1,1)=squeeze(derivatives2{j}(:,1,1,2,2))
    %D_uv F^j(0,0)=squeeze(derivatives2{j}(:,1,2,1,1))
    %D_uv F^j(0,1)=squeeze(derivatives2{j}(:,1,2,1,2))
    %D_uv F^j(1,0)=squeeze(derivatives2{j}(:,1,2,2,1))
    %D_uv F^j(1,1)=squeeze(derivatives2{j}(:,1,2,2,2))  
    %D_vv F^j(0,0)=squeeze(derivatives2{j}(:,2,2,1,1))
    %D_vv F^j(0,1)=squeeze(derivatives2{j}(:,2,2,1,2))
    %D_vv F^j(1,0)=squeeze(derivatives2{j}(:,2,2,2,1))
    %D_vv F^j(1,1)=squeeze(derivatives2{j}(:,2,2,2,2))
end

%Construction of CC_edges
for iref = 1:numel(interfaces_all)
    clear patch
    if ~isempty(interfaces_all(iref).patch1)
        patch(1) = interfaces_all(iref).patch1; %LEFT
        side(1) = interfaces_all(iref).side1; %LEFT
    end
    if ~isempty(interfaces_all(iref).patch2)
        patch(2) = interfaces_all(iref).patch2; %RIGHT
        side(2) = interfaces_all(iref).side2; %RIGHT
    end
    nnz_el=find(patch>0);
    if side(1)==3 && side(2)==1
        patch=flip(patch);
        side=flip(side);
    end
  
  %STEP 2 - Stuff necessary to evaluate the geo_mapping and its derivatives
  for ii = nnz_el % The two patches (L-R)
  %in this cycle we must add to the already existing CC_vertices, the part of the "discarded" interface functions    
    brk = cell (1,msh.ndim);
    knots = space.sp_patch{patch(ii)}.knots;
    %Greville points for G^1 conditions system
    geo_knot_v1=geometry(patch(ii)).nurbs.knots{1};
    p1=geometry(patch(ii)).nurbs.order(1);
    aug_geo_knot1=[geo_knot_v1(1)*ones(1,p1+1)  repelem(geo_knot_v1(p1+1:end-p1-1),3)  geo_knot_v1(end)*ones(1,p1+1)];
    geo_knot_v2=geometry(patch(ii)).nurbs.knots{2};
    p2=geometry(patch(ii)).nurbs.order(2);
    aug_geo_knot2=[geo_knot_v2(1)*ones(1,p2+1)  repelem(geo_knot_v2(p2+1:end-p2-1),3)  geo_knot_v2(end)*ones(1,p2+1)];
    grev_pts{1}=aveknt(aug_geo_knot1,p1+1);
    grev_pts{2}=aveknt(aug_geo_knot2,p2+1);
    for idim = 1:msh.ndim
        if (numel(grev_pts{idim}) > 1)
            brk{idim} = [knots{idim}(1), grev_pts{idim}(1:end-1) + diff(grev_pts{idim})/2, knots{idim}(end)];
        else
            brk{idim} = [knots{idim}(1) knots{idim}(end)];
        end
    end  
    msh_grev = msh_cartesian (brk, grev_pts, [], geometry(patch(ii)), 'boundary', true, 'der2',false);
    msh_side_int{ii} = msh_boundary_side_from_interior (msh_grev, side(ii));
    msh_side_int{ii} = msh_precompute (msh_side_int{ii});      
    geo_map_jac{ii} = msh_side_int{ii}.geo_map_jac; %rdim x ndim x 1 x n_grev_pts (rdim->physical space, ndim->parametric space)
%     if ii==1 && (side(ii)==3 || side(ii)==4)
%         disp('done1')
%         geo_map_jac{ii}=flip(geo_map_jac{ii});
%     elseif ii==2 && (side(ii)==1 || side(ii)==2)
%         disp('done2')
%         geo_map_jac{ii}=flip(geo_map_jac{ii});
%     end
  end
  
  if length(patch)==2 %as it is placed now, this check computes the matrices corresponding only to interior edges
  %STEP 3 - Assembling and solving G^1 conditions system  %this must depend on orientation!
  if side(2)==1 || side(2)==2
      v=grev_pts{2}(:);
  else
      v=grev_pts{1}(:);
  end
  ngrev=numel(v);
  DuFR_x=reshape(geo_map_jac{1}(1,1,:,:),ngrev,1); %column vector
  DuFR_y=reshape(geo_map_jac{1}(2,1,:,:),ngrev,1); %column vector
  DvFL_x=reshape(geo_map_jac{2}(1,2,:,:),ngrev,1); %column vector
  DvFL_y=reshape(geo_map_jac{2}(2,2,:,:),ngrev,1); %column vector
  DvFR_x=reshape(geo_map_jac{1}(1,2,:,:),ngrev,1); %column vector
  DvFR_y=reshape(geo_map_jac{1}(2,2,:,:),ngrev,1); %column vector
  
  A_full=[(1-v).*DvFL_x v.*DvFL_x (1-v).*DuFR_x v.*DuFR_x (1-v).^2.*DvFR_x 2*(1-v).*v.*DvFR_x v.^2.*DvFR_x;...
     (1-v).*DvFL_y v.*DvFL_y (1-v).*DuFR_y v.*DuFR_y (1-v).^2.*DvFR_y 2*(1-v).*v.*DvFR_y v.^2.*DvFR_y];
 if rank(A_full)==6
     A=A_full(:,2:end);
     b=-A_full(:,1);
     sols=A\b;
     alpha0_n(2)=1; %R
     alpha1_n(2)=sols(1); %R
     alpha0_n(1)=sols(2); %L
     alpha1_n(1)=sols(3); %L
     beta0_n=sols(4);
     beta1_n=sols(5);
     beta2_n=sols(6);
 else
     A=A_full(:,3:end);
     b=-sum(A_full(:,1:2),2);
     sols=A\b;
     alpha0_n(2)=1; %R
     alpha1_n(2)=1; %R
     alpha0_n(1)=sols(1); %L
     alpha1_n(1)=sols(2); %L
     beta0_n=sols(3);
     beta1_n=sols(4);
     beta2_n=sols(5);     
 end
 
 %keyboard
 %STEP 4 - Normalizing the alphas
 %C1=((alpha1_n(1)-alpha0_n(1))^2)/3+((alpha1_n(2)-alpha0_n(2))^2)/3 + (alpha1_n(1)-alpha0_n(1))*alpha0_n(1)+(alpha1_n(2)-alpha0_n(2))*alpha0_n(2)...
 %   +alpha0_n(1)^2+alpha0_n(2)^2;
 %C2=(alpha1_n(1)-alpha0_n(1))-(alpha1_n(2)-alpha0_n(2))+2*alpha0_n(1)-2*alpha0_n(2);
 %gamma=-C2/(2*C1);
 C1=alpha0_n(1)^2+alpha0_n(1)*alpha1_n(1)+alpha1_n(1)^2+alpha0_n(2)^2+alpha0_n(2)*alpha1_n(2)+alpha1_n(2)^2;
 C2=alpha0_n(1)+alpha1_n(1)+alpha0_n(2)+alpha1_n(2);
 gamma=3*C2/(2*C1);
 alpha0(2)=alpha0_n(2)*gamma; %R
 alpha1(2)=alpha1_n(2)*gamma; %R
 alpha0(1)=alpha0_n(1)*gamma; %L
 alpha1(1)=alpha1_n(1)*gamma; %L
 bbeta0=beta0_n*gamma;
 bbeta1=beta1_n*gamma;
 bbeta2=beta2_n*gamma;
 

 %STEP 5 - Computing the betas
 %alphas and beta evaluated at 0,1,1/2
 alpha_L_0=alpha0(1); %alpha_L(0)
 alpha_L_1=alpha1(1); %alpha_L(1)
 alpha_L_12=(alpha0(1)+alpha1(1))/2; %alpha_L(1/2)
 alpha_R_0=alpha0(2); %alpha_L(0)
 alpha_R_1=alpha1(2); %alpha_L(1)
 alpha_R_12=(alpha0(2)+alpha1(2))/2; %alpha_L(1/2)  
 beta_0=bbeta0; %beta(0)
 beta_1=bbeta2; %beta(1)
 beta_12=(bbeta0+bbeta2)/4+bbeta1/2; %beta(1/2)
 
  %Computing the matrix of the system considering the relationship between beta^L, beta^R and beta
 M=[alpha_L_0 0 -alpha_R_0 0; 0 alpha_L_1 0 -alpha_R_1; alpha_L_12/2 alpha_L_12/2 -alpha_R_12/2 -alpha_R_12/2];
 
 if rank(M)==3
     
 %Computing beta1_L, beta0_R, beta1_R in terms of beta0_L
 quant1=(alpha_R_12/2-(alpha_R_0*alpha_L_12)/(2*alpha_L_0))/((alpha_R_1*alpha_L_12)/(2*alpha_L_1)-alpha_R_12/2);
 quant2=(beta_12-(beta_0*alpha_L_12)/(2*alpha_L_0)-(beta_1*alpha_L_12)/(2*alpha_L_1))/((alpha_R_1*alpha_L_12)/(2*alpha_L_1)-alpha_R_12/2); 
 
 %beta1_L=a+b*beta0_L,  beta0_R=c+d*beta0_L,  beta1_R=e+f*beta0_L, where
 a=quant2; b=quant1;
 c=beta_0/alpha_L_0; d=alpha_R_0/alpha_L_0;
 e=(beta_1+alpha_R_1*quant2)/alpha_L_1; f=alpha_R_1*quant1/alpha_L_1;
 
 %We determine beta0_L by minimizing the sum of the norms of beta_L and beta_R
 C1=((b-1)^2)/3+(b-1)+((f-d)^2)/3+(f-d)*d+d^2+1;
 C2=2*a*(b-1)/3+a+2*(e-c)*(f-d)/3+(e-c)*d+(f-d)*c+2*c*d;
 beta0(1)=-C2/(2*C1); %L
 beta1(1)=a+b*beta0(1); %L
 beta0(2)=c+d*beta0(1); %R
 beta1(2)=e+f*beta0(1); %R
 
 else
     
 %Computing beta0_R in terms of beta0_L and beta1_R in terms of beta1_L: 
 %beta0_R=a+b*beta0_L,  beta1_R=c+d*beta1_L, where
 a=beta_0/alpha_L_0; b=alpha_R_0/alpha_L_0;
 c=beta_1/alpha_L_1; d=alpha_R_1/alpha_L_1;
 
 %We determine beta0_L and beta_1_L by minimizing the sum of the norms of beta_L and beta_R
 %The resuting system is
 M2=[2*(1+b^2) 1+b*d; 1+b*d 2*(1+d^2)];
 M2b=[-b*c-2*a*b; -a*d-2*c*d];
 sol=M2\M2b;
 beta0(1)= sol(1); %L
 beta1(1)= sol(2); %L
 beta0(2)= a+b*beta0(1); %R
 beta1(2)= c+d*beta1(1); %R
 
 end
 
 %Saving alphas and betas (first column=L, second column=R)
 all_alpha0(iref,:)=alpha0;
 all_alpha1(iref,:)=alpha1;
 all_beta0(iref,:)=beta0;
 all_beta1(iref,:)=beta1;  
    
% Compute the Greville points, and the auxiliary mesh and space objects
    for ii = 1:2 % The two patches (R-L)
      brk = cell (1,msh.ndim);
      knots = space.sp_patch{patch(ii)}.knots;
      
% The knot vectors for the N0 and N1 basis functions, as in Mario's notation
% Only univariate knot vectors are computed
%%    %ind1  = [2 2 1 1]; ind2 = [1 1 2 2]
      ind2 = ceil (side(ii)/2);
      ind1 = setdiff (1:msh.ndim, ind2);

      degree = space.sp_patch{patch(1)}.degree(ind1);
      knots0 = space.knots0_patches{patch(ii)}{ind1};
      knots1 = space.knots1_patches{patch(ii)}{ind1};
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

% Degree and first length in the direction normal to the interface
      ind2 = ceil (side(ii)/2);
      degu = space.sp_patch{patch(ii)}.degree(ind2);
      knt = unique (space.sp_patch{patch(ii)}.knots{ind2});
      if (mod (side(ii), 2) == 1)
        tau1 = knt(2) - knt(1);
      else
        tau1 = knt(end) - knt(end-1);
      end

% For now we assume that the orientation is as in the paper, and we do not need any reordering
%XXXX    [sp_bnd(2), msh_side(2)] = reorder_elements_and_quad_points (sp_bnd(2), msh_side(2), interfaces(iref), ndim);
% msh_side contains only the univariate parametrization of the boundary (dependence on u)
% msh_side_int contains information for the bivariate parametrization (dependence on u and v)
% sp_grev contains the value and derivatives of the basis functions, at the Greville points
      msh_side(ii) = msh_eval_boundary_side (msh_grev, side(ii));
      msh_side_int{ii} = msh_boundary_side_from_interior (msh_grev, side(ii));
      sp_aux = space.sp_patch{patch(ii)}.constructor (msh_side_int{ii});
      msh_side_int{ii} = msh_precompute (msh_side_int{ii});
%       sp_grev = sp_precompute_param (sp_aux, msh_side_int{ii}, 'value', true, 'gradient', true); % This could be done in just one element

% The univariate spaces for the basis functions N^{p,r+1} (knots0) and N^{p-1,r} (knots1)
      sp0 = sp_bspline (knots0, degree, msh_grev.boundary(side(ii)));
      sp1 = sp_bspline (knots1, degree-1, msh_grev.boundary(side(ii)));
      sp0_struct = sp_precompute_param (sp0, msh_grev.boundary(side(ii)), 'value', true, 'gradient', true);
      sp1_struct = sp_precompute_param (sp1, msh_grev.boundary(side(ii)), 'value', true, 'gradient', true);

% The univariate space for the basis functions N^{p,r} on the interface
%       spn = space.sp_patch{patch(ii)}.boundary(side(ii)).constructor(msh_grev.boundary(side(ii)))
      ind = [2 2 1 1];
      knotsn = knots{ind(side(ii))};
      spn = sp_bspline (knotsn, degree, msh_grev.boundary(side(ii)));
      spn_struct = sp_precompute_param (spn, msh_grev.boundary(side(ii)), 'value', true, 'gradient', true);

% Matrix for the linear systems, (14)-(16) in Mario's notes
      A = sparse (msh_side(ii).nel, msh_side(ii).nel);
      for jj = 1:msh_side(ii).nel
        A(jj,spn_struct.connectivity(:,jj)) = spn_struct.shape_functions(:,:,jj);
      end

%alphas and betas
        alpha{ii}=abs(alpha0(ii)*(1-grev_pts{2}')+alpha1(ii)*grev_pts{2}'); %we have to take the absolute value to make it work for any orientation
        beta{ii}=beta0(ii)*(1-grev_pts{2}')+beta1(ii)*grev_pts{2}';       

% RHS for the first linear system, (14) in Mario's notes
      rhss = sparse (msh_side(ii).nel, sp0_struct.ndof);
      for jj = 1:msh_side(ii).nel
        rhss(jj,sp0_struct.connectivity(:,jj)) = sp0_struct.shape_functions(:,:,jj);
      end
      coeff0{ii} = A \ rhss;
      coeff0{ii}(abs(coeff0{ii}) < 1e-12) = 0; % Make more sparse

% RHS for the second linear system, (15) in Mario's notes
      rhsb = sparse (msh_side(ii).nel, sp0_struct.ndof);
      if (side(ii) == 1) %paper case
        val_grad = sp_aux.sp_univ(1).shape_function_gradients(2);
      elseif (side(ii) == 2)
        val_grad = sp_aux.sp_univ(1).shape_function_gradients(end-1);
      elseif (side(ii) == 3)
        val_grad = sp_aux.sp_univ(2).shape_function_gradients(2);
      elseif (side(ii) == 4)
        val_grad = sp_aux.sp_univ(2).shape_function_gradients(end-1);
      end
      val = val_grad * (tau1 / degu)^2;
      for jj = 1:msh_side(ii).nel  %paper case - check beta
        val_aux = -val * beta{ii}(jj); %added a minus, as in (6) of multipatch paper
        rhsb(jj,sp0_struct.connectivity(:,jj)) = sp0_struct.shape_function_gradients(:,:,:,jj) * val_aux;
      end
      rhsb = rhsb + rhss;
      coeff1{ii} = A \ rhsb;
      coeff1{ii}(abs(coeff1{ii}) < 1e-12) = 0; % Make more sparse
      
% RHS for the third linear system, (16) in Mario's notes
% We need this change of sign to make it work for general orientation.
% I don't understand why we don't need it in the previous system.
      val_grad = val_grad * (-1)^(side(ii)+1);

      rhsc = sparse (msh_side(ii).nel, sp1_struct.ndof);
      val = val_grad* (tau1 / degu)^2;  %WARNING: WE DIVIDED BY tau1/degu, which REQUIRES A SMALL MODIFICATION IN THE REF MASK (ADD MULT. BY 1/2) 
      for jj = 1:msh_side(ii).nel
        val_aux = val * alpha{ii}(jj)* (-1)^(ii-1); %with the multipatch settings must be multiplied by -1 for left patch;
        rhsc(jj,sp1_struct.connectivity(:,jj)) = sp1_struct.shape_functions(:,:,jj) * val_aux;
      end
      coeff2{ii} = A \ rhsc;
      coeff2{ii}(abs(coeff2{ii}) < 1e-12) = 0; % Make more sparse
      
% Pass the coefficients to the tensor product basis
% The numbering (ndof) only works for the two patch case, for now
      ndof = sp0_struct.ndof + sp1_struct.ndof;
      CC_edges{ii,iref} = sparse (space.ndof_per_patch(patch(ii)), ndof);
      
      ndof_dir = space.sp_patch{patch(ii)}.ndof_dir;
      if (side(ii) == 1)
        ind0 = sub2ind (ndof_dir, ones(1,spn.ndof), 1:spn.ndof);
        ind1 = sub2ind (ndof_dir, 2*ones(1,spn.ndof), 1:spn.ndof);
      elseif (side(ii) == 2)
        ind0 = sub2ind (ndof_dir, ndof_dir(1) * ones(1,spn.ndof), 1:spn.ndof);
        ind1 = sub2ind (ndof_dir, (ndof_dir(1)-1) * ones(1,spn.ndof), 1:spn.ndof);
      elseif (side(ii) == 3)
        ind0 = sub2ind (ndof_dir, 1:spn.ndof, ones(1,spn.ndof));
        ind1 = sub2ind (ndof_dir, 1:spn.ndof, 2*ones(1,spn.ndof));
      elseif (side(ii) == 4)
        ind0 = sub2ind (ndof_dir, 1:spn.ndof, ndof_dir(2) * ones(1,spn.ndof));
        ind1 = sub2ind (ndof_dir, 1:spn.ndof, (ndof_dir(2)-1) * ones(1,spn.ndof));
      end
      if (ii == 2 & interfaces_all(iref).ornt == -1) %this should be still the same for the multipatch
        ind0 = fliplr (ind0);
        ind1 = fliplr (ind1);
        ind0_s{iref}=ind0; %saving the indices for later use;
        ind1_s{iref}=ind1;
        coeff0{2} = flipud (fliplr (coeff0{2}));
        coeff1{2} = flipud (fliplr (coeff1{2}));
        coeff2{2} = flipud (fliplr (coeff2{2}));
      end

      CC_edges{ii,iref}(ind0,1:sp0_struct.ndof) = coeff0{ii};  %coefficients of trace basis functions...
      CC_edges{ii,iref}(ind1,1:sp0_struct.ndof) = coeff1{ii};
      CC_edges{ii,iref}(ind1,sp0_struct.ndof+1:ndof) = coeff2{ii};%CORRECT REMOVING THIS? * (-1)^(ii+1); %...and derivative basis functions (associated to iref-th interface)
      
      %keeping the part of the "actually active" edge functions, the remaining part saved in CC_edges_discarded
      CC_edges_discarded{ii,iref}=CC_edges{ii,iref}(:,[1 2 3 sp0_struct.ndof-2:sp0_struct.ndof+2 ndof-1 ndof]); %dimension: n^2 x 10
      CC_edges{ii,iref}=CC_edges{ii,iref}(:,[4:sp0_struct.ndof-3 sp0_struct.ndof+3:ndof-2]);
    end
  else 
     n0=n-k; 
     ndof=n0+n-k-1;
     CC_edges{ii,iref} = sparse (space.ndof_per_patch(patch(ii)), ndof);
%      size(CC_edges{ii,iref})
%      keyboard
     CC_edges_discarded{ii,iref}=CC_edges{ii,iref}(:,[1 2 3 n0-3:n0+1 ndof-1 ndof]); %dimension: n^2 x 10
     CC_edges{ii,iref}=CC_edges{ii,iref}(:,[4:n0-4 n0+2:ndof-2]);
  end
    
%CHECKING G^1 condition  
% bbeta=alpha{1}.*beta{2}-alpha{2}.*beta{1};
% geo_map_jac{1} = msh_side_int{1}.geo_map_jac; %rdim x ndim x 1 x n_grev_pts (rdim->physical space, ndim->parametric space)
% geo_map_jac{2} = msh_side_int{2}.geo_map_jac; %rdim x ndim x 1 x n_grev_pts (rdim->physical space, ndim->parametric space)
%   for ii = 1:2 % The two patches (L-R)
%     brk = cell (1,msh.ndim);
%     knots = space.sp_patch{patch(ii)}.knots;
%     for idim = 1:msh.ndim
%         %Greville points for G^1 conditions system
%         geo_knot_v1=geometry(patch(ii)).nurbs.knots{1};
%         p1=geometry(patch(ii)).nurbs.order(1);
%         aug_geo_knot1=[geo_knot_v1(1)*ones(1,p1+1)  repelem(geo_knot_v1(p1+1:end-p1-1),3)  geo_knot_v1(end)*ones(1,p1+1)];
%         geo_knot_v2=geometry(patch(ii)).nurbs.knots{2};
%         p2=geometry(patch(ii)).nurbs.order(2);
%         aug_geo_knot2=[geo_knot_v2(1)*ones(1,p2+1)  repelem(geo_knot_v2(p2+1:end-p2-1),3)  geo_knot_v2(end)*ones(1,p2+1)];
%         grev_pts{1}=aveknt(aug_geo_knot1,p1+1);
%         grev_pts{2}=aveknt(aug_geo_knot2,p2+1);
%         if (numel(grev_pts{idim}) > 1)
%             brk{idim} = [knots{idim}(1), grev_pts{idim}(1:end-1) + diff(grev_pts{idim})/2, knots{idim}(end)];
%         else
%             brk{idim} = [knots{idim}(1) knots{idim}(end)];
%         end
%     end  
%   
%     msh_grev = msh_cartesian (brk, grev_pts, [], geometry(patch(ii)), 'boundary', true, 'der2',false);
%     msh_side_int{ii} = msh_boundary_side_from_interior (msh_grev, side(ii));
%     msh_side_int{ii} = msh_precompute (msh_side_int{ii});      
%     geo_map_jac{ii} = msh_side_int{ii}.geo_map_jac; %rdim x ndim x 1 x n_grev_pts (rdim->physical space, ndim->parametric space) 
%   end
% alpha{1}=(alpha0(1)*(1-grev_pts{2}')+alpha1(1)*grev_pts{2}'); %we have to take the absolute value to make it work for any orientation
% beta{1}=beta0(1)*(1-grev_pts{2}')+beta1(1)*grev_pts{2}';
% alpha{2}=(alpha0(2)*(1-grev_pts{2}')+alpha1(2)*grev_pts{2}'); %we have to take the absolute value to make it work for any orientation
% beta{2}=beta0(2)*(1-grev_pts{2}')+beta1(2)*grev_pts{2}';  
% bbeta=alpha{1}.*beta{2}-alpha{2}.*beta{1};
% DuFL_x=squeeze(geo_map_jac{1}(1,1,:,:)); %column vector
% DuFL_y=squeeze(geo_map_jac{1}(2,1,:,:)); %column vector
% DuFR_x=squeeze(geo_map_jac{2}(1,1,:,:)); %column vector
% DuFR_y=squeeze(geo_map_jac{2}(2,1,:,:)); %column vector
% DvFL_x=squeeze(geo_map_jac{1}(1,2,:,:)); %column vector
% DvFL_y=squeeze(geo_map_jac{1}(2,2,:,:)); %column vector
%alpha{1}.*beta{2}-alpha{2}.*beta{1}-bbeta
% alpha{2}.*DuFL_x-alpha{1}.*DuFR_x+bbeta.*DvFL_x
% alpha{2}.*DuFL_y-alpha{1}.*DuFR_y+bbeta.*DvFL_y
%   A_full=[(1-v).*DuFL_x v.*DuFL_x -(1-v).*DuFR_x -v.*DuFR_x (1-v).^2.*DvFL_x 2*(1-v).*v.*DvFL_x v.^2.*DvFL_x;...
%      (1-v).*DuFL_y v.*DuFL_y -(1-v).*DuFR_y -v.*DuFR_y (1-v).^2.*DvFL_y 2*(1-v).*v.*DvFL_y v.^2.*DvFL_y];
% alpha{1}
% alpha{2}
% beta{1}
% beta{2}
% bbeta
% alpha0
% alpha1
% beta0
% beta1
%pause
% grev_pts{1}
% grev_pts{2}

end



%We assume that the local numbering of interfaces and patches is such that
%vertices(kver).interface(im) is the interface between
%vertices(kver).patches(im) and vertices(kver).patches(im+1)
MM=cell(2,numel(vertices));
V=cell(numel(vertices),1);
E=cell(numel(vertices),1);
%Previously:
%MM=cell(2,nu,numel(vertices));
%V=cell(nu,numel(vertices));
%E=cell(nu+1,numel(vertices));
sides=[1 4;2 3;1 2;4 3]; %on i-th row the indices of the endpoints of the i-th side (bottom-up, left-right)

for kver=1:numel(vertices)
    
    %Everything must be updated by using interfaces_all instead of interfaces TO DO
    ver_patches=[];%vector with indices of patches containing the vertex
    ver_patches_nabla={}; %cell array containing jacobians
    ver_ind=[];%vector containing local index of vertex in the patch
    nu=numel(vertices(kver).interfaces);
    
    boundary_v=0;
    for im=1:nu
        inter=vertices(kver).interfaces(im);
        if isempty(interfaces_all(inter).side1) | isempty(interfaces_all(inter).side2)
            boundary_v=1;
            break
        end
    end
    if boundary_v==0
    
%     for h=1:numel(vertices(kver).interfaces)
%         hint=vertices(kver).interfaces(h);
%         ver_patches=[ver_patches interfaces(hint).patch1 interfaces(hint).patch2];    
%         ver_ind=[ver_ind sides(interfaces(hint).side1,vertices(kver).ind)...
%                  sides(interfaces(hint).side2,vertices(kver).ind)];
%     end
%     ver_patches=unique(ver_patches,'stable');
%     ver_ind=unique(ver_ind,'stable');
    
    for im=1:nu %cycle over all the interfaces containing the vertex
        inter=vertices(kver).interfaces(im); %global index of the interface 
        patch_ind1=interfaces_all(inter).patch1; %global index of left patch of im-th interface
        patch_ind2=interfaces_all(inter).patch2; %global index of right patch of im-th interface
        vertex_ind1=sides(interfaces_all(inter).side1,vertices(kver).ind(im)); %local index of vertex in left patch
        vertex_ind2=sides(interfaces_all(inter).side2,vertices(kver).ind(im)); %local index of vertex in right patch
        ver_patches=[ver_patches patch_ind1 patch_ind2];
        ver_ind=[ver_ind vertex_ind1 vertex_ind2];
        %compute t(0) and t'(0), d(0) and d'(0)
        switch vertex_ind1
            case 1 %vertex (0,0)
                Du_F00=squeeze(derivatives1{patch_ind1}(:,1,1));
                Dv_F00=squeeze(derivatives1{patch_ind1}(:,2,1));
                Duv_F00=squeeze(derivatives2{patch_ind1}(:,1,2,1));
                Dvv_F00=squeeze(derivatives2{patch_ind1}(:,2,2,1));
            case 2 %vertex (0,1)
                Du_F00=squeeze(derivatives1{patch_ind1}(:,1,2));
                Dv_F00=-squeeze(derivatives1{patch_ind1}(:,2,2));
                Duv_F00=-squeeze(derivatives2{patch_ind1}(:,1,2,2));
                Dvv_F00=squeeze(derivatives2{patch_ind1}(:,2,2,2));  
            case 3 %vertex (1,0)
                Du_F00=-squeeze(derivatives1{patch_ind1}(:,1,3));
                Dv_F00=squeeze(derivatives1{patch_ind1}(:,2,3));
                Duv_F00=-squeeze(derivatives2{patch_ind1}(:,1,2,3));
                Dvv_F00=squeeze(derivatives2{patch_ind1}(:,2,2,3));  
            case 4 %vertex (1,1)
                Du_F00=-squeeze(derivatives1{patch_ind1}(:,1,4));
                Dv_F00=-squeeze(derivatives1{patch_ind1}(:,2,4));
                Duv_F00=squeeze(derivatives2{patch_ind1}(:,1,2,4));
                Dvv_F00=squeeze(derivatives2{patch_ind1}(:,2,2,4));
        end
        %Store the jacobian of F for the left patch
        ver_patches_nabla{2*im-1}=[Du_F00 Dv_F00];
        
        t0(im,:)=Dv_F00;
        t0p(im,:)=Dvv_F00;
        d0(im,:)=(Du_F00+(all_beta0(inter,1)*(1-0)+all_beta1(inter,1)*0)*Dv_F00)...
            /(all_alpha0(inter,1)*(1-0)+all_alpha1(inter,1)*0);
        d0p(im,:)=(-all_alpha0(inter,1)*(Du_F00+(all_beta0(inter,1)*(1-0)+all_beta1(inter,1)*0)*Dv_F00)+...
                   (all_alpha0(inter,1)*(1-0)+all_alpha1(inter,1)*0)*(Duv_F00...
                   -all_beta0(inter,1)*Dv_F00+(all_beta0(inter,1)*(1-0)+all_beta1(inter,1)*0)*Dvv_F00...
                   ))/(all_alpha0(inter,1)*(1-0)+all_alpha1(inter,1)*0)^2;  
        mix_der2(im,:)=Duv_F00;
        %We need to get the jacobian also for the right patch
        switch vertex_ind2
            case 1 %vertex (0,0)
                Du_F00=squeeze(derivatives1{patch_ind2}(:,1,1));
                Dv_F00=squeeze(derivatives1{patch_ind2}(:,2,1));
            case 2 %vertex (0,1)
                Du_F00=squeeze(derivatives1{patch_ind2}(:,1,2));
                Dv_F00=-squeeze(derivatives1{patch_ind2}(:,2,2));
            case 3 %vertex (1,0)
                Du_F00=-squeeze(derivatives1{patch_ind2}(:,1,3));
                Dv_F00=squeeze(derivatives1{patch_ind2}(:,2,3));
            case 4 %vertex (1,1)
                Du_F00=-squeeze(derivatives1{patch_ind2}(:,1,4));
                Dv_F00=-squeeze(derivatives1{patch_ind2}(:,2,4));
        end
        ver_patches_nabla{2*im}=[Du_F00 Dv_F00];
        
        %Pick the correct part of CC_edges_discarded %TO BE FIXED
        if vertices(kver).ind(im)==1
            E{kver}{im,1}=CC_edges_discarded{1,inter}(:,[1 2 3 7 8]); %part of the matrix corresponding to edge functions close to the vertex
            E{kver}{im,2}=CC_edges_discarded{2,inter}(:,[1 2 3 7 8]);
        else
            E{kver}{im,1}=CC_edges_discarded{1,inter}(:,[4 5 6 9 10]);
            E{kver}{im,2}=CC_edges_discarded{2,inter}(:,[4 5 6 9 10]);
        end    
    end
    [ver_patches, ind_patch_sigma, ind_patch_rep]=unique(ver_patches,'stable');
    %ver_ind=unique(ver_ind,'rows','stable');
    
    %if the number of patches coincides with the number of interfaces, 
    %we add one fictional interface coinciding with the first one
    %(just for coding-numbering reasons)
    if numel(ver_patches)==nu
        t0(nu+1,:)=t0(1,:);
        t0p(nu+1,:)=t0p(1,:);
        d0(nu+1,:)=d0(1,:);
        d0p(nu+1,:)=d0p(1,:);
        E{kver}{nu,1}=E{kver}{1,1};
        E{kver}{nu,2}=E{kver}{1,2};
    end
    
    %computing sigma
    sigma=0;
    for im=1:nu
        sigma=sigma+norm(ver_patches_nabla{ind_patch_sigma(im)},Inf);
    end
    sigma=1/(sigma/(p*(k+1)*nu));
    
    %computing matrices MM and V
    for im=1:nu
        %assemble matrix (not final: Ms and Vs, then updated with the "discarded parts" of edge functions)
        n1=space.sp_patch{ver_patches(im)}.ndof_dir(1); %dimension of tensor-product space in the patch (dir 1)
        n2=space.sp_patch{ver_patches(im)}.ndof_dir(2); %dimension of tensor-product space in the patch (dir 2)
        im_edges=ceil(find(ind_patch_rep==im)/2); %global indices of the edges 
        im1=im_edges(1); im2=im_edges(2);
        j=1;
        for j1=0:2
            for j2=0:2-j1 %the following computations work in the standard case
                d00=(j1==0)*(j2==0);
                %M_{i_{m-1},i}
                d10_a=[(j1==1)*(j2==0), (j1==0)*(j2==1)]*t0(im1,:)';
                d20_a=t0(im1,:)*[(j1==2)*(j2==0), (j1==1)*(j2==1); (j1==1)*(j2==1), (j1==0)*(j2==2)]*t0(im1,:)'+...
                      [(j1==1)*(j2==0), (j1==0)*(j2==1)]*t0p(im1,:)';
                d01_a=[(j1==1)*(j2==0), (j1==0)*(j2==1)]*d0(im1,:)';
                d11_a=t0(im1,:)*[(j1==2)*(j2==0), (j1==1)*(j2==1); (j1==1)*(j2==1), (j1==0)*(j2==2)]*d0(im1,:)'+...
                      [(j1==1)*(j2==0), (j1==0)*(j2==1)]*d0p(im1,:)';
                MM{1,kver}{im}(:,j)=sigma^(j1+j2)*[d00, d00+d10_a/(p*(k+1)), d00+2*d10_a/(p*(k+1))+d20_a/(p*(p-1)*(k+1)^2),...
                                                 d01_a/(p*(k+1)), d01_a/(p*(k+1))+d11_a/(p*(p-1)*(k+1)^2)]';
                %M_{i_{m+1},i}
                d10_b=[(j1==1)*(j2==0), (j1==0)*(j2==1)]*t0(im2,:)';
                d20_b=t0(im2,:)*[(j1==2)*(j2==0), (j1==1)*(j2==1); (j1==1)*(j2==1), (j1==0)*(j2==2)]*t0(im2,:)'+...
                      [(j1==1)*(j2==0), (j1==0)*(j2==1)]*t0p(im2,:)';
                d01_b=[(j1==1)*(j2==0), (j1==0)*(j2==1)]*d0(im2,:)';
                d11_b=t0(im2,:)*[(j1==2)*(j2==0), (j1==1)*(j2==1); (j1==1)*(j2==1), (j1==0)*(j2==2)]*d0(im2,:)'+...
                      [(j1==1)*(j2==0), (j1==0)*(j2==1)]*d0p(im2,:)';                
                MM{2,kver}{im}(:,j)=sigma^(j1+j2)*[d00, d00+d10_b/(p*(k+1)), d00+2*d10_b/(p*(k+1))+d20_b/(p*(p-1)*(k+1)^2),...
                                                 d01_b/(p*(k+1)), d01_b/(p*(k+1))+d11_b/(p*(p-1)*(k+1)^2)]';     
                %V_{i_m,i}  
                d11_c=t0(im1,:)*[(j1==2)*(j2==0), (j1==1)*(j2==1); (j1==1)*(j2==1), (j1==0)*(j2==2)]*t0(im2,:)'+...
                      [(j1==1)*(j2==0), (j1==0)*(j2==1)]*mix_der2(im1,:)';
                V{kver}{im}=zeros(n1*n2,6);
                V{kver}{im}([1, 2, n2+1, n2+2],j)=[d00, d00+d10_b/(p*(k+1)), d00+d10_a/(p*(k+1)),...
                                                  d00+ (d10_a+d10_b+d11_c/(p*(k+1)))/(p*(k+1))]';    
                j=j+1;
            end
        end
        if interfaces_all(vertices(kver).interfaces(im1)).side2==ver_patches(im)%the considered patch is to the right of edge im1
            E1=E{kver}{im1,2};
        else
            E1=E{kver}{im1,1};
        end
        if interfaces_all(vertices(kver).interfaces(im2)).side2==ver_patches(im)%the considered patch is to the right of edge im2
            E2=E{kver}{im2,2};
        else
            E2=E{kver}{im2,1};
        end
        CC_vertices{ver_patches(im),kver} = E1*MM{1,kver}{im} + E2*MM{2,kver}{im} - V{kver}{im};%E{kver}{im1,2}*MM{1,kver}{im} + E{kver}{im2,1}*MM{2,kver}{im} - V{kver}{im};
    end
    end

end

end

%TO DO:
%- case of boundary edges

%TO BE TESTED
%- number of functions
%- plots of the functions
%- continuity 