% SP_MULTIPATCH_C1: Constructor of the class for multipatch spaces with C1 continuity
%
%     sp = sp_multipatch_C1 (spaces, msh, geometry, edges, vertices)
%
% INPUT:
%
%    spaces:   cell-array of space objects, one for each patch (see sp_scalar)
%    msh:      mesh object that defines the multipatch mesh (see msh_multipatch)
%    geometry: geometry struct (see mp_geo_load)
%    edges:    struct array similar to interfaces (see vertices_struct)
%    vertices: struct array with vertex information (see vertices_struct)
%
% OUTPUT:
%
%    sp: object representing the discrete function space of vector-valued functions, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                          DESCRIPTION
%        npatch          (scalar)                        number of patches
%        ncomp           (scalar)                        number of components of the functions of the space (always equal to one)
%        ndof            (scalar)                        total number of degrees of freedom after gluing patches together
%        ndof_per_patch  (1 x npatch array)              number of degrees of freedom per patch, without gluing
%        ndof_interior_per_patch  (1 x npatch array)     number of interior basis functions on each patch
%        ndof_interior   (scalar)                        total number of interior degrees of freedom
%        ndof_edges      (scalar)                        total number of edge degrees of freedom
%        ndof_vertices   (scalar)                        total number of vertex degrees of freedom
%        dofs_on_edge    (1 x nedges cell-array)         global numbering of edge functions on each edge
%        dofs_on_vertex  (1 x nvertices cell-array)      global numbering of the six functions on each vertex
%        sp_patch        (1 x npatch cell-array)         the input spaces, one space object for each patch (see sp_scalar)
%        CC_edges        (2 x nedges cell-array)         matrices for B-spline representation of edge basis functions
%        CC_vertices     (1 x nvertices cell-array)      matrices for B-spline representation of vertex basis functions, 
%                                                         as many as the valence for each vertex
%        constructor     function handle                 function handle to construct the same discrete space in a different msh
%
%       METHODS
%       Methods that give a structure with all the functions computed in a certain subset of the mesh
%         sp_evaluate_element_list: compute basis functions (and derivatives) in a given list of elements
%
%       Methods for post-processing, that require a computed vector of degrees of freedom
%         sp_l2_error:    compute the error in L2 norm
%         sp_h1_error:    compute the error in H1 norm
%         sp_h2_error:    compute the error in H2 norm (only for planar surfaces)
%         sp_eval_phys:   compute the value of the solution in a given set of points
%         sp_to_vtk:      export the computed solution to a pvd file, using a Cartesian grid of points on each patch
%         sp_plot_solution: plot the solution in Matlab (only for the scalar-valued case)
%
%       Methods for basic connectivity operations
%         sp_get_basis_functions: compute the functions that do not vanish in a given list of elements
%         sp_get_cells:           compute the cells on which a list of functions do not vanish
%         sp_get_neighbors:       compute the neighbors, functions that share at least one element with a given one
%         sp_get_functions_on_patch: compute the indices of non-vanishing C^1 functions on a patch
%         sp_get_local_interior_functions: compute the local indices of interior functions on a patch
%         sp_get_vertex_neighbors: compute the indices of elements in the support of vertex functions
%
%       Other methods
%         sp_refine: generate a refined space, and subdivision matrices for the univariate spaces
%         sp_compute_Cpatch: compute the matrix for B-splines representation on a patch
%         sp_compute_Cpatch_vector: compute the matrix for B-splines representation on a patch
%                                   for vector-valued spaces with the same space on each component
%
% Copyright (C) 2015-2024 Rafael Vazquez
% Copyright (C) 2019-2023 Cesare Bracco
% Copyright (C) 2022 Andrea Farahat
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

function sp = sp_multipatch_C1 (spaces, msh, geometry, interfaces_all, vertices)

  if (~all (cellfun (@(x) isa (x, 'sp_scalar'), spaces)))
    error ('All the spaces in the array should be of the sp_scalar class')
  end
  ndof_all = cellfun (@(x) x.ndof, spaces);
  degree_all = cellfun (@(x) x.degree, spaces, 'UniformOutput', false);
  degree_all = cell2mat (degree_all.');

  sp.npatch = numel (spaces);
  if (sp.npatch ~= msh.npatch)
    error ('The list of spaces does not correspond to the mesh')
  end

  if (~all (all (degree_all == spaces{1}.degree, 2), 1))
    error ('All the patches should have the same degree, at least for now.')
  end
  
  for iptc = 1:numel(geometry)
    knots = spaces{iptc}.knots;
    breaks = cellfun (@unique, knots, 'UniformOutput', false);
    mult = cellfun (@histc, knots, breaks, 'UniformOutput', false);
    if (any ([mult{:}] < 2))
      error ('The regularity should be at most degree minus two')
    end
    for idim = 1:msh.ndim
      if (any (mult{idim}(2:end-1) > spaces{iptc}.degree(idim) - 1))
        error ('The regularity should not be lower than one')
      end
    end
    if (numel(breaks{1}) ~= numel(breaks{2}))
      error ('The number of internal knots should be the same in every patch and every coordinate, at least for now.')
    end
  end

  sp.ncomp = spaces{1}.ncomp;
  sp.transform = spaces{1}.transform;
  
  for iptc = 1:sp.npatch
    if (~strcmpi (spaces{iptc}.transform, 'grad-preserving'))
      error ('The transform to the physical domain should be the same for all the spaces, and the grad-preserving one')
    end
    if (~strcmpi (spaces{iptc}.space_type, 'spline'))
      error ('C1 continuity is only implemented for splines, not for NURBS')
    end
  end

  sp.ndof = 0;
  sp.ndof_per_patch = ndof_all;
  sp.sp_patch = spaces;
  
% Knot vectors of auxiliary spaces
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

% Computation of the number of degrees of freedom, and the coefficients to
% express them as linear combinations of B-splines.
% See the function compute_coefficients below for details.  
  sp.ndof_interior = 0;
  for iptc = 1:sp.npatch
    [ii,jj] = ndgrid (3:spaces{iptc}.ndof_dir(1)-2, 3:spaces{iptc}.ndof_dir(2)-2);
    interior_dofs = sub2ind (spaces{iptc}.ndof_dir, ii(:)', jj(:)');
    
    sp.ndof_interior = sp.ndof_interior + numel (interior_dofs);
    sp.ndof_interior_per_patch(iptc) = numel (interior_dofs);
  end
  clear interior_dofs

  warnaux = warning('query','geopdes:nrbmultipatch');
  warning ('off', 'geopdes:nrbmultipatch')
  [ndof_per_interface, CC_edges, ndof_per_vertex, CC_vertices, v_fun_matrices] = ...
    compute_coefficients (sp, msh, geometry, interfaces_all, vertices);
  warning (warnaux.state, 'geopdes:nrbmultipatch')

%Sigmas, K and V matrices used in the construction of vertex functions
  sp.vertex_function_matrices = v_fun_matrices;

% Total number of functions, and of functions of each type
  sp.ndof_edges = sum (ndof_per_interface);
  sp.ndof_vertices = sum (ndof_per_vertex);
  sp.ndof = sp.ndof_interior + sp.ndof_edges + sp.ndof_vertices;
  
  sp.dofs_on_edge = cell (1, numel(interfaces_all));
  for intrfc = 1:numel(interfaces_all)
    global_indices = sp.ndof_interior + sum(ndof_per_interface(1:intrfc-1)) + (1:ndof_per_interface(intrfc));
    sp.dofs_on_edge{intrfc} = global_indices;
  end

  sp.dofs_on_vertex = cell (1, numel(vertices));
  for ivrt = 1:numel(vertices)
    global_indices = sp.ndof_interior + sp.ndof_edges + sum(ndof_per_vertex(1:ivrt-1)) + (1:ndof_per_vertex(ivrt));
    sp.dofs_on_vertex{ivrt} = global_indices;
  end
  
  sp.CC_edges = CC_edges;
  sp.CC_vertices = CC_vertices;
  
% Store information about the geometry, for simplicity
  sp.interfaces = interfaces_all;
  sp.vertices = vertices;
  sp.geometry = geometry;
  
  sp.constructor = @(MSH) sp_multipatch_C1 (patches_constructor(spaces, MSH), MSH, geometry, interfaces_all, vertices);
    function spaux = patches_constructor (spaces, MSH)
      for ipatch = 1:MSH.npatch
        spaux{ipatch} = spaces{ipatch}.constructor(MSH.msh_patch{ipatch});
      end
    end

  sp = class (sp, 'sp_multipatch_C1');  
end


function [ndof_per_interface, CC_edges, ndof_per_vertex, CC_vertices, v_fun_matrices] = compute_coefficients (space, msh, geometry, interfaces_all, vertices)
% OUTPUT from compute_coefficients  
% ndof_per_interface: number of edge functions on each interface, array of size 1 x numel(interfaces);
% CC_edges: cell array of size 2 x numel(interfaces), the two corresponds
%   to the two patches on the interface. The matrix CC_edges{ii,jj} has size
%      sp.ndof_per_patch(patch) x ndof_per_interface(jj)
%      with patch = interfaces(jj).patches(ii);
% ndof_per_vertex: number of vertex functions on each vertex. An array of size numel(vertices). Even if the value is always six.
% CC_vertices: cell array of size 1 x numel(vertices), each entry containin another 
%   cell-array of size 1 x vertices(ii).valence_p (local number of patches)
%   The matrix CC_vertices{jj}{ii} has size sp.ndof_per_patch(patch) x ndof_per_vertex(ii)
%      with patch = vertices(jj).patches(ii).

% Initialize cell array containing sigma, K and V matrices for each vertex
v_fun_matrices = cell(2, numel(vertices));

%Initialize output variables with correct size
ndof_per_interface = zeros (1, numel(interfaces_all));
ndof_per_vertex = 6*ones(1,numel(vertices));

CC_edges = cell (2, numel (interfaces_all));
CC_edges_discarded = cell (2, numel (interfaces_all));
CC_vertices = cell (1, numel(vertices));
for jj = 1:numel(vertices)
  CC_vertices{jj} = cell (1,vertices(jj).valence_p);
  for kk=1:vertices(jj).valence_p
    CC_vertices{jj}{kk}=zeros(space.sp_patch{kk}.ndof,6);
  end
end

deg_aux = space.sp_patch{1}.degree(1);
breaks_m = unique (space.sp_patch{1}.knots{1});
nbrk_aux = numel(breaks_m) - 1; % This is k+1 in the paper
mult = histc (space.sp_patch{1}.knots{1},breaks_m);
% reg_aux = deg_aux - max(mult(2:end-1));
reg_aux = deg_aux - max(mult([2 end-1]));

% Check that regularity of the first knot is the same in all patches and directions
for iptc = 1:space.npatch
  for idim = 1:msh.ndim
    breaks_m = cellfun (@unique, space.sp_patch{iptc}.knots, 'UniformOutput', false);
    mult = histc(space.sp_patch{iptc}.knots{1},breaks_m{1});
    reg{idim} = space.sp_patch{iptc}.degree(1) - max(mult([2 end-1]));
  end
  if (any (reg{1} ~= reg{2}))
    error ('The regularity of the first internal knot should always be the same, at least for now.')
  end
end

all_alpha0 = zeros(numel(interfaces_all),2);
all_alpha1 = zeros(numel(interfaces_all),2);
all_beta0 = zeros(numel(interfaces_all),2);
all_beta1 = zeros(numel(interfaces_all),2);

%Computation of CC_edges
for iref = 1:numel(interfaces_all)
  patches = [interfaces_all(iref).patch1 interfaces_all(iref).patch2];
  npatches_on_edge = numel (patches);
  operations = interfaces_all(iref).operations;
% Auxiliary geometry with orientation as in the paper
  geo_local = reorientation_patches (operations, geometry(patches));
  sides = [1 3];

% Compute gluing data  
  geo_map_jac = cell (npatches_on_edge, 1);
  for ii = 1:npatches_on_edge
    brk = cell (1,msh.ndim); 
    grev_pts = cell (1, msh.ndim);
    knt = geo_local(ii).nurbs.knots;
    order = geo_local(ii).nurbs.order;
    for idim = 1:msh.ndim
      grev_pts{idim} = [knt{idim}(order(idim)) (knt{idim}(order(idim))+knt{idim}(end-order(idim)+1))/2 knt{idim}(end-order(idim)+1)];
      brk{idim} = [knt{idim}(order(idim)), grev_pts{idim}(1:end-1) + diff(grev_pts{idim})/2, knt{idim}(end-order(idim)+1)];
    end
    msh_grev = msh_cartesian (brk, grev_pts, [], geo_local(ii), 'boundary', true, 'der2',false);
    msh_side_interior = msh_boundary_side_from_interior (msh_grev, sides(ii));
    msh_side_interior = msh_precompute (msh_side_interior);
    geo_map_jac{ii} = msh_side_interior.geo_map_jac; %rdim x ndim x 1 x n_grev_pts (rdim->physical space, ndim->parametric space)
  end
    
  if (msh.ndim + 1 == msh.rdim)
    [alpha0, alpha1, beta0, beta1] = compute_gluing_data_surf (geo_map_jac, grev_pts, sides);
  else
    [alpha0, alpha1, beta0, beta1] = compute_gluing_data (geo_map_jac, grev_pts, sides);
  end
  clear geo_map_jac msh_grev msh_side_interior grev_pts

 %Saving alphas and betas (first column=i_0, second column=i_1)
  all_alpha0(iref,:) = alpha0;
  all_alpha1(iref,:) = alpha1;
  all_beta0(iref,:) = beta0;
  all_beta1(iref,:) = beta1;  

% Compute the Greville points, and the auxiliary mesh and space objects for
%  functions with reduced degree or increased regularity
  for ii = 1:npatches_on_edge
%    %ind1  = [2 2 1 1]; ind2 = [1 1 2 2]
    ind2 = ceil (sides(ii)/2);
    ind1 = setdiff (1:msh.ndim, ind2);

    brk = cell (1,msh.ndim);
    grev_pts = cell (1, msh.ndim);
    degrees = degree_reorientation (space.sp_patch{patches(ii)}.degree, operations(ii,3));
    degree = degrees(ind1);

    knots = knot_vector_reorientation (space.sp_patch{patches(ii)}.knots, operations(ii,:));
    knots0 = knot_vector_reorientation (space.knots0_patches{patches(ii)}, operations(ii,:));
    knots0 = knots0{ind1};
    knots1 = knot_vector_reorientation (space.knots1_patches{patches(ii)}, operations(ii,:));
    knots1 = knots1{ind1};
    for idim = 1:msh.ndim
      grev_pts{idim} = aveknt (knots{idim}, degrees(idim)+1);%space.sp_patch{patch(ii)}.degree(idim)+1); 
      grev_pts{idim} = grev_pts{idim}(:)';
      brk{idim} = [knots{idim}(1), grev_pts{idim}(1:end-1) + diff(grev_pts{idim})/2, knots{idim}(end)];
    end
    msh_grev = msh_cartesian (brk, grev_pts, [], geo_local(ii), 'boundary', true, 'der2',false);

% Degree and first length in the direction normal to the interface
    degu = degrees(ind2);
    knt = unique (knots{ind2});
    tau1 = knt(2) - knt(1);

% sp_aux contains the value and derivatives of the basis functions, at the Greville points
    msh_side = msh_eval_boundary_side (msh_grev, sides(ii));
    msh_side_interior = msh_boundary_side_from_interior (msh_grev, sides(ii));
    sp_aux = sp_bspline (knots, degrees, msh_side_interior);

% Univariate spaces for the basis functions N^{p,r+1} (knots0) and N^{p-1,r} (knots1) and N^{p,r} on the interface
    sp0 = sp_bspline (knots0, degree, msh_grev.boundary(sides(ii)));
    sp1 = sp_bspline (knots1, degree-1, msh_grev.boundary(sides(ii)));
    sp0_struct = sp_precompute_param (sp0, msh_grev.boundary(sides(ii)), 'value', true, 'gradient', true);
    sp1_struct = sp_precompute_param (sp1, msh_grev.boundary(sides(ii)), 'value', true, 'gradient', true);

    knotsn = knots{ind1};
    spn = sp_bspline (knotsn, degree, msh_grev.boundary(sides(ii)));
    spn_struct = sp_precompute_param (spn, msh_grev.boundary(sides(ii)), 'value', true, 'gradient', true);

% Matrix for the linear systems
    A = sparse (msh_side.nel, msh_side.nel);
    for jj = 1:msh_side.nel
      A(jj,spn_struct.connectivity(:,jj)) = spn_struct.shape_functions(:,:,jj);
    end

% alphas and betas (gluing data)
    alpha = alpha0(ii)*(1-grev_pts{3-ii}') + alpha1(ii)*grev_pts{3-ii}';
    beta = beta0(ii)*(1-grev_pts{3-ii}') + beta1(ii)*grev_pts{3-ii}';

% RHS and solution of the linear systems, (14)-(16) in Mario's notes
    rhss = sparse (msh_side.nel, sp0_struct.ndof);
    for jj = 1:msh_side.nel
      rhss(jj,sp0_struct.connectivity(:,jj)) = sp0_struct.shape_functions(:,:,jj);
    end
    coeff0 = A \ rhss;
    coeff0(abs(coeff0) < 1e-12) = 0; % Make more sparse

    rhsb = sparse (msh_side.nel, sp0_struct.ndof);
    if (ii == 1) % paper case, (sides(ii) == 1)
      val_grad = sp_aux.sp_univ(1).shape_function_gradients(2);
    elseif (ii == 2) %The other paper case, (sides(ii) == 3)
      val_grad = sp_aux.sp_univ(2).shape_function_gradients(2);
    end
    val = val_grad * (tau1 / degu)^2;
    for jj = 1:msh_side.nel
      val_aux = -val * beta(jj);
      rhsb(jj,sp0_struct.connectivity(:,jj)) = sp0_struct.shape_function_gradients(:,:,:,jj) * val_aux;
    end
    rhsb = rhsb + rhss;
    coeff1 = A \ rhsb;
    coeff1(abs(coeff1) < 1e-12) = 0; % Make more sparse

    rhsc = sparse (msh_side.nel, sp1_struct.ndof);
    val = val_grad* (tau1 / degu);
    for jj = 1:msh_side.nel
      val_aux = val * alpha(jj)* (-1)^(ii+1);
      rhsc(jj,sp1_struct.connectivity(:,jj)) = sp1_struct.shape_functions(:,:,jj) * val_aux;
    end
    coeff2 = A \ rhsc;
    coeff2(abs(coeff2) < 1e-12) = 0; % Make more sparse

% Pass the coefficients to the tensor product basis
    ndof_dir = sp_aux.ndof_dir; %space.sp_patch{patch(ii)}.ndof_dir;
    ndof_dir_original = space.sp_patch{patches(ii)}.ndof_dir;
    if (ii == 1)     %(sides(ii) == 1)
      ind0 = sub2ind (ndof_dir, ones(1,spn.ndof), 1:spn.ndof);
      ind1 = sub2ind (ndof_dir, 2*ones(1,spn.ndof), 1:spn.ndof);
    elseif (ii == 2) %(sides(ii) == 3)
      ind0 = sub2ind (ndof_dir, 1:spn.ndof, ones(1,spn.ndof));
      ind1 = sub2ind (ndof_dir, 1:spn.ndof, 2*ones(1,spn.ndof));
    end
    indices_reoriented = indices_reorientation (ndof_dir_original, operations(ii,:));
    ind0_ornt = indices_reoriented(ind0);
    ind1_ornt = indices_reoriented(ind1);

% Store the coefficients in CC_edges, and compute the number of functions
    ndof_edge = sp0_struct.ndof + sp1_struct.ndof;
    trace_functions = 4:sp0_struct.ndof-3;
    deriv_functions = sp0_struct.ndof + (3:sp1_struct.ndof-2);
    active_functions = union (trace_functions, deriv_functions);
    discarded_functions = setdiff (1:ndof_edge, active_functions);

    CC = sparse (space.ndof_per_patch(patches(ii)), ndof_edge);
    CC(ind0_ornt,1:sp0_struct.ndof) = coeff0;
    CC(ind1_ornt,1:sp0_struct.ndof) = coeff1;
    CC(ind1_ornt,sp0_struct.ndof+(1:sp1_struct.ndof)) = coeff2;

    ndof_per_interface(iref) = numel(active_functions);    
    CC_edges{ii,iref} = sparse (space.ndof_per_patch(patches(ii)), ndof_per_interface(iref));
    CC_edges{ii,iref} = CC(:,active_functions);
    CC_edges_discarded{ii,iref} = CC(:,discarded_functions);
  end
end

% Computation of CC_vertices
% Auxiliary constants to deal with different regularity cases in the M matrices
if (reg_aux < deg_aux-2)
  const1 = 2; const2 = 1;
else
  const1 = 3; const2 = 2;
end

for kver = 1:numel(vertices)
  edges = vertices(kver).edges;
  patches = vertices(kver).patches;
  valence_e = vertices(kver).valence_e;
  valence_p = vertices(kver).valence_p;
  operations = vertices(kver).patch_reorientation;
  edge_orientation = vertices(kver).edge_orientation;
  
  %Initializing the cell-array containing, for vertex kver, K and V matrices
  KV_matrices=cell(1,valence_p);
  
  geo_local = reorientation_patches (operations, geometry(patches));
  % Here we need to further modify geo_local by rotating the parametrization
  %  of the patches, that is, for each ip in patches, rotate the control points
  %  geometry(ip).nurbs.coefs (4 x (number of control points in dir 1) x (number of control points in dir 2) matrix)

% Precompute the derivatives and compute sigma
  sigma = 0;
  for iptc = 1:valence_p
    knots = space.sp_patch{iptc}.knots;
    for idim = 1:msh.ndim
      brk{idim} = [knots{idim}(1) knots{idim}(end)];
    end
    msh_pts_der1 = msh_cartesian (brk, {0 0}, [], geo_local(iptc),'boundary', true, 'der2', true);
    msh_der = msh_precompute (msh_pts_der1);
    derivatives_new1{iptc} = msh_der.geo_map_jac; %rdim x ndim x (n_pts{1}x n_pts{2}) (rdim->physical space, ndim->parametric space)
    derivatives_new2{iptc} = msh_der.geo_map_der2; %rdim x ndim x ndim x n_pts{1} x n_pts{2}
    
    sigma = sigma + norm (derivatives_new1{iptc}, Inf);
  end
  sigma = deg_aux * nbrk_aux * valence_p / sigma;
% Store sigma for output
  v_fun_matrices{1,kver} = sigma;
  
  if (msh.ndim+1 == msh.rdim)
    % Tangent vectors
    Du_F = derivatives_new1{1}(:,1);
    Dv_F = derivatives_new1{1}(:,2);
    % Normal vector
    normal = cross(Du_F,Dv_F);
    unit_normal = normal/norm(normal);
    % Vector orthogonal to n and z along which the geometry is rotated (rotation axis)
    r = cross(unit_normal,[0 0 1]');
    % Angle between n and z
    cos_th = [0 0 1]*unit_normal;
    sin_th = norm(r);
    % Rotation matrix
    R = [cos_th + r(1)^2*(1-cos_th),r(1)*r(2)*(1-cos_th)-r(3)*sin_th,r(1)*r(3)*(1-cos_th)+r(2)*sin_th;...
         r(1)*r(2)*(1-cos_th)+r(3)*sin_th,cos_th + r(2)^2*(1-cos_th),r(2)*r(3)*(1-cos_th)-r(1)*sin_th;...
         r(1)*r(3)*(1-cos_th)-r(2)*sin_th,r(3)*r(2)*(1-cos_th)+r(1)*sin_th,cos_th + r(3)^2*(1-cos_th)];
     
    for ipatch = 1:valence_p
      coeff_hom = geo_local(ipatch).nurbs.coefs;
      coeff_w = coeff_hom(4,:,:);
      coeff_eu = coeff_hom(1:3,:,:)./(repmat(coeff_w,3,1,1)); %Euclidean coefficients
      rot_coeff = ones(4,size(coeff_eu,2),size(coeff_eu,3));
      for k = 1:size(coeff_eu,3)
        rot_coeff(1:3,:,k) = R*(coeff_eu(:,:,k)-coeff_eu(:,1,1))+coeff_eu(:,1,1);
      end
      rot_nrb(ipatch) = geometry(patches(ipatch)).nurbs;
      rot_nrb(ipatch).coefs = rot_coeff;
    end
   
    rot_geo = mp_geo_load(rot_nrb);
      
    for iptc = 1:valence_p
      knots = space.sp_patch{iptc}.knots;
      for idim = 1:msh.ndim
        brk{idim} = [knots{idim}(1) knots{idim}(end)];
      end
      rot_msh_pts_der1 = msh_cartesian (brk, {0 0}, [], rot_geo(iptc),'boundary', true, 'der2', true);
      rot_msh_der = msh_precompute (rot_msh_pts_der1);
      derivatives_new1{iptc} = rot_msh_der.geo_map_jac; %rdim x ndim x (n_pts{1}x n_pts{2}) (rdim->physical space, ndim->parametric space)
      derivatives_new2{iptc} = rot_msh_der.geo_map_der2; %rdim x ndim x ndim x n_pts{1} x n_pts{2}
    end
  end
  
  for ipatch = 1:valence_p
    prev_edge = ipatch;
    next_edge = mod(ipatch, valence_e) + 1;

% Compute gluing data, and edge functions from CC_edges_discarded
    if (edge_orientation(prev_edge) == 1)
      alpha_prev = all_alpha0(edges(prev_edge),2);
      beta_prev = all_beta0(edges(prev_edge),2);
      alpha_der_prev = -all_alpha0(edges(prev_edge),2) + all_alpha1(edges(prev_edge),2);
      beta_der_prev = -all_beta0(edges(prev_edge),2) + all_beta1(edges(prev_edge),2);
      E_prev = CC_edges_discarded{2,edges(prev_edge)}(:,[1 2 3 7 8]);
    else
      alpha_prev = -all_alpha1(edges(prev_edge),1);
      beta_prev = -all_beta1(edges(prev_edge),1);
      alpha_der_prev = -all_alpha0(edges(prev_edge),1) + all_alpha1(edges(prev_edge),1);
      beta_der_prev = -all_beta0(edges(prev_edge),1) + all_beta1(edges(prev_edge),1);
      E_prev = CC_edges_discarded{1,edges(prev_edge)}(:,[6 5 4 10 9]);
    end
    if (edge_orientation(next_edge) == 1)
      alpha_next = all_alpha0(edges(next_edge),1);
      beta_next = all_beta0(edges(next_edge),1);
      alpha_der_next = -all_alpha0(edges(next_edge),1) + all_alpha1(edges(next_edge),1);
      beta_der_next = -all_beta0(edges(next_edge),1) + all_beta1(edges(next_edge),1);
      E_next = CC_edges_discarded{1,edges(next_edge)}(:,[1 2 3 7 8]);
    else
      alpha_next = -all_alpha1(edges(next_edge),2);
      beta_next = -all_beta1(edges(next_edge),2);
      alpha_der_next = -all_alpha0(edges(next_edge),2) + all_alpha1(edges(next_edge),2);
      beta_der_next = -all_beta0(edges(next_edge),2) + all_beta1(edges(next_edge),2);
      E_next = CC_edges_discarded{2,edges(next_edge)}(:,[6 5 4 10 9]);
    end
    
    Du_F = derivatives_new1{ipatch}(1:2,1);
    Dv_F = derivatives_new1{ipatch}(1:2,2);
    Duu_F = derivatives_new2{ipatch}(1:2,1,1);
    Duv_F = derivatives_new2{ipatch}(1:2,1,2);
    Dvv_F = derivatives_new2{ipatch}(1:2,2,2);
    
% Edge information
    t0_prev = Du_F;
    t0_next = Dv_F;
    t0p_prev = Duu_F;
    t0p_next = Dvv_F;

    d0_prev = -(Dv_F + beta_prev * Du_F) / alpha_prev;
    d0_next =  (Du_F + beta_next * Dv_F) / alpha_next;
    d0p_prev = ( alpha_der_prev * (Dv_F + beta_prev*Du_F) - ...
                 alpha_prev * (Duv_F + beta_der_prev*Du_F + beta_prev*Duu_F)) / alpha_prev^2;
    d0p_next = (-alpha_der_next * (Du_F + beta_next*Dv_F) + ...
                 alpha_next * (Duv_F + beta_der_next*Dv_F + beta_next*Dvv_F)) / alpha_next^2;

% Compute M and V matrices
    ndof = space.sp_patch{patches(ipatch)}.ndof;
    K_prev = sparse (5,6); K_next = sparse (5,6);
    VV = sparse (ndof,6);
    
    ndof_dir = space.sp_patch{patches(ipatch)}.ndof_dir;
    all_indices = indices_reorientation (ndof_dir, operations(ipatch,:));
    corner_4dofs = all_indices(1:2,1:2);
    jfun = 1;
    for j1 = 0:2
      for j2 = 0:2-j1
        mat_deltas = [(j1==2)*(j2==0), (j1==1)*(j2==1); (j1==1)*(j2==1), (j1==0)*(j2==2)];
        vec_deltas = [(j1==1)*(j2==0); (j1==0)*(j2==1)];
        c0 = (j1==0)*(j2==0);
        c1_a = vec_deltas.'*t0_prev;
        c2_a = t0_prev.'*mat_deltas*t0_prev + vec_deltas.'*t0p_prev;
        d0_a = vec_deltas.'*d0_prev;
        d1_a = t0_prev.'*mat_deltas*d0_prev + vec_deltas.'*d0p_prev;
 
        c1_b = vec_deltas.'*t0_next;
        c2_b = t0_next.'*mat_deltas*t0_next + vec_deltas.'*t0p_next;
        d0_b = vec_deltas.'*d0_next;
        d1_b = t0_next.'*mat_deltas*d0_next + vec_deltas.'*d0p_next;

        K_prev(:,jfun) = sigma^(j1+j2)*[c0, ...
                                        c0+c1_a/(deg_aux*nbrk_aux), ...
                                        c0+const1*c1_a/(deg_aux*nbrk_aux)+const2*c2_a/(deg_aux*(deg_aux-1)*nbrk_aux^2), ...
                                        d0_a/(deg_aux*nbrk_aux), ...
                                        d0_a/(deg_aux*nbrk_aux)+d1_a/(deg_aux*(deg_aux-1)*nbrk_aux^2)].';
        K_next(:,jfun) = sigma^(j1+j2)*[c0, ...
                                        c0+c1_b/(deg_aux*nbrk_aux), ...
                                        c0+const1*c1_b/(deg_aux*nbrk_aux)+const2*c2_b/(deg_aux*(deg_aux-1)*nbrk_aux^2), ...
                                        d0_b/(deg_aux*nbrk_aux), ...
                                        d0_b/(deg_aux*nbrk_aux)+d1_b/(deg_aux*(deg_aux-1)*nbrk_aux^2)].';
        e11 = t0_prev.'*mat_deltas*t0_next + vec_deltas.'*Duv_F;
        VV(corner_4dofs,jfun) = sigma^(j1+j2)*[c0, ...
                                               c0+c1_a/(deg_aux*nbrk_aux), ...
                                               c0+c1_b/(deg_aux*nbrk_aux), ...
                                               c0+(c1_a+c1_b+e11/(deg_aux*nbrk_aux))/(deg_aux*nbrk_aux)]'; 
        jfun = jfun+1;
      end
    end
    CC_vertices{kver}{ipatch} = E_prev*K_prev + E_next*K_next - VV;
% Store the matrices for output
    KV_matrices{ipatch}.K_prev = K_prev;
    KV_matrices{ipatch}.K_next = K_next;
    KV_matrices{ipatch}.V = VV;
    KV_matrices{ipatch}.E_prev = E_prev;
    KV_matrices{ipatch}.E_next = E_next;
  end
  v_fun_matrices{2,kver}=KV_matrices;
end


end

function [alpha0, alpha1, beta0, beta1] = compute_gluing_data (geo_map_jac, grev_pts, side)
% OUTPUT from gluing data
% The gluing data are computed as alpha = alpha0*(1-t) + alpha1*t.
% The two components of alpha0, alpha1, beta0, beta1 refer to the neighboring
%  patches: first or second ((i,0) and (i,1) in the draft, L and R in old notation).

% For boundary edges we set alpha=1 and beta=0. For simplicity, we store an inexistent second patch.
  if (numel (geo_map_jac) == 1)
    alpha0 = [1 1];
    alpha1 = [1 1];
    beta0 = [0 0]; 
    beta1 = [0 0];
    return
  end

  % Assemble and solve G^1 conditions system
  if (side(2)==1 || side(2)==2)
    v = grev_pts{2}(:);
  else
    v = grev_pts{1}(:);
  end
  ngrev = numel(v);
  DuFR_x = reshape (geo_map_jac{1}(1,1,:,:), ngrev, 1);
  DuFR_y = reshape (geo_map_jac{1}(2,1,:,:), ngrev, 1);
  DvFL_x = reshape (geo_map_jac{2}(1,2,:,:), ngrev, 1);
  DvFL_y = reshape (geo_map_jac{2}(2,2,:,:), ngrev, 1);
  DvFR_x = reshape (geo_map_jac{1}(1,2,:,:), ngrev, 1);
  DvFR_y = reshape (geo_map_jac{1}(2,2,:,:), ngrev, 1);
  
  A_full = [(1-v).*DvFL_x v.*DvFL_x (1-v).*DuFR_x v.*DuFR_x (1-v).^2.*DvFR_x 2*(1-v).*v.*DvFR_x v.^2.*DvFR_x;...
       (1-v).*DvFL_y v.*DvFL_y (1-v).*DuFR_y v.*DuFR_y (1-v).^2.*DvFR_y 2*(1-v).*v.*DvFR_y v.^2.*DvFR_y];
  if (rank(A_full)==6)
    A = A_full(:,2:end);
    b = -A_full(:,1);
    sols = A\b;
    alpha0_n(1) = 1;
    alpha1_n(1) = sols(1);
    alpha0_n(2) = sols(2);
    alpha1_n(2) = sols(3);
    beta0_n = sols(4);
    beta1_n = sols(5);
    beta2_n = sols(6);
  else
    A = A_full(:,3:end); % FIX: not a square matrix
    b = -sum(A_full(:,1:2),2);
    sols = A\b;
    alpha0_n(1) = 1;
    alpha1_n(1) = 1;
    alpha0_n(2) = sols(1);
    alpha1_n(2) = sols(2);
    beta0_n = sols(3);
    beta1_n = sols(4);
    beta2_n = sols(5);     
  end
 
 % Normalize the alphas
  C1 = alpha0_n(2)^2+alpha0_n(2)*alpha1_n(2)+alpha1_n(2)^2+alpha0_n(1)^2+alpha0_n(1)*alpha1_n(1)+alpha1_n(1)^2;
  C2 = alpha0_n(2)+alpha1_n(2)+alpha0_n(1)+alpha1_n(1);
  gamma = 3*C2/(2*C1);
  alpha0(1) = alpha0_n(1)*gamma;
  alpha1(1) = alpha1_n(1)*gamma;
  alpha0(2) = alpha0_n(2)*gamma;
  alpha1(2) = alpha1_n(2)*gamma;
  bbeta0 = beta0_n*gamma;
  bbeta1 = beta1_n*gamma;
  bbeta2 = beta2_n*gamma;
 
 % Compute the betas, evaluating alphas and betas at 0, 1/2, 1.
  alpha_R_0 = alpha0(1);
  alpha_R_1 = alpha1(1);
  alpha_R_12 = (alpha0(1)+alpha1(1))/2;
  alpha_L_0 = alpha0(2);
  alpha_L_1 = alpha1(2);
  alpha_L_12 = (alpha0(2)+alpha1(2))/2;
  beta_0 = bbeta0;
  beta_1 = bbeta2;
  beta_12 = (bbeta0+bbeta2)/4+bbeta1/2;
 
 % Compute the matrix of the system considering the relationship between beta^(i,0), beta^(i,1) and beta
  M = [alpha_R_0 0 alpha_L_0 0; ...
       0 alpha_R_1 0 alpha_L_1; ...
       alpha_R_12/2 alpha_R_12/2 alpha_L_12/2 alpha_L_12/2];
 
  if (rank(M)==3)
 % Compute beta1_R, beta0_L, beta1_L in terms of beta0_R
    quant1 = (-alpha_L_12/2 + (alpha_L_0*alpha_R_12)/(2*alpha_R_0)) / ...
      (-(alpha_L_1*alpha_R_12)/(2*alpha_R_1) + alpha_L_12/2);
    quant2 = (beta_12-(beta_0*alpha_R_12)/(2*alpha_R_0) - (beta_1*alpha_R_12)/(2*alpha_R_1)) / ...
      (-(alpha_L_1*alpha_R_12)/(2*alpha_R_1) + alpha_L_12/2); 
 
 % beta1_R=a+b*beta0_R,  beta0_L=c+d*beta0_R,  beta1_L=e+f*beta0_R, where
    a = quant2; b = quant1;
    c = beta_0/alpha_R_0; d = -alpha_L_0/alpha_R_0;
    e = (beta_1-alpha_L_1*quant2)/alpha_R_1; f = -alpha_L_1*quant1/alpha_R_1;
      
 %We determine beta0_R by minimizing the sum of the norms of beta_R and beta_L
    C1 = ((b-1)^2)/3 + (b-1) + ((f-d)^2)/3 + (f-d)*d + d^2 + 1;
    C2 = 2*a*(b-1)/3 + a + 2*(e-c)*(f-d)/3 + (e-c)*d + (f-d)*c + 2*c*d;
    beta0(1) = -C2/(2*C1);
    beta1(1) = a + b*beta0(1);
    beta0(2) = c + d*beta0(1);
    beta1(2) = e + f*beta0(1);

  else
 % Compute beta0_L in terms of beta0_R and beta1_L in terms of beta1_R:
 % beta0_R=a+b*beta0_L,  beta1_R=c+d*beta1_L, where
    a = beta_0/alpha_R_0; b = -alpha_L_0/alpha_R_0;
    c = beta_1/alpha_R_1; d = -alpha_L_1/alpha_R_1;
 
 % Determine beta0_R and beta_1_R by minimizing the sum of the norms of beta_R and beta_L
    M2 = [2*(1+b^2) 1+b*d; 1+b*d 2*(1+d^2)];
    M2b = [-b*c-2*a*b; -a*d-2*c*d];
    sol = M2\M2b;
    beta0(1)= sol(1);
    beta1(1)= sol(2);
    beta0(2)= a + b*beta0(1);
    beta1(2)= c + d*beta1(1);
  end 
  
end


function [alpha0, alpha1, beta0, beta1] = compute_gluing_data_surf (geo_map_jac, grev_pts, side)
% OUTPUT from gluing data
% The gluing data are computed as alpha = alpha0*(1-t) + alpha1*t.
% The two components of alpha0, alpha1, beta0, beta1 refer to the neighboring
%  patches: first or second ((i,0) and (i,1) in the draft, L and R in old notation).

% For boundary edges we set alpha=1 and beta=0. For simplicity, we store an inexistent second patch.
  if (numel (geo_map_jac) == 1)
    alpha0 = [1 1];
    alpha1 = [1 1];
    beta0 = [0 0]; 
    beta1 = [0 0];
    return
  end

  % Assemble and solve G^1 conditions system
  if (side(2)==1 || side(2)==2)
    v = grev_pts{2}(:);
  else
    v = grev_pts{1}(:);
  end
  ngrev = numel(v);
  DuFR_x = reshape (geo_map_jac{1}(1,1,:,:), ngrev, 1);
  DuFR_y = reshape (geo_map_jac{1}(2,1,:,:), ngrev, 1);
  DuFR_z = reshape (geo_map_jac{1}(3,1,:,:), ngrev, 1);
  DvFL_x = reshape (geo_map_jac{2}(1,2,:,:), ngrev, 1);
  DvFL_y = reshape (geo_map_jac{2}(2,2,:,:), ngrev, 1);
  DvFL_z = reshape (geo_map_jac{2}(3,2,:,:), ngrev, 1);
  DvFR_x = reshape (geo_map_jac{1}(1,2,:,:), ngrev, 1);
  DvFR_y = reshape (geo_map_jac{1}(2,2,:,:), ngrev, 1);
  DvFR_z = reshape (geo_map_jac{1}(3,2,:,:), ngrev, 1);
  
  A_full = [(1-v).*DvFL_x v.*DvFL_x (1-v).*DuFR_x v.*DuFR_x (1-v).^2.*DvFR_x 2*(1-v).*v.*DvFR_x v.^2.*DvFR_x;...
       (1-v).*DvFL_y v.*DvFL_y (1-v).*DuFR_y v.*DuFR_y (1-v).^2.*DvFR_y 2*(1-v).*v.*DvFR_y v.^2.*DvFR_y;
       (1-v).*DvFL_z v.*DvFL_z (1-v).*DuFR_z v.*DuFR_z (1-v).^2.*DvFR_z 2*(1-v).*v.*DvFR_z v.^2.*DvFR_z];
   
  if (rank(A_full)==6)
    A = A_full(:,2:end);
    b = -A_full(:,1);
    sols = A\b;
    alpha0_n(1) = 1;
    alpha1_n(1) = sols(1);
    alpha0_n(2) = sols(2);
    alpha1_n(2) = sols(3);
    beta0_n = sols(4);
    beta1_n = sols(5);
    beta2_n = sols(6);
  else
    A = A_full(:,3:end); % FIX: not a square matrix
    b = -sum(A_full(:,1:2),2);
    sols = A\b;
    alpha0_n(1) = 1;
    alpha1_n(1) = 1;
    alpha0_n(2) = sols(1);
    alpha1_n(2) = sols(2);
    beta0_n = sols(3);
    beta1_n = sols(4);
    beta2_n = sols(5);     
  end
 
 % Normalize the alphas
  C1 = alpha0_n(2)^2+alpha0_n(2)*alpha1_n(2)+alpha1_n(2)^2+alpha0_n(1)^2+alpha0_n(1)*alpha1_n(1)+alpha1_n(1)^2;
  C2 = alpha0_n(2)+alpha1_n(2)+alpha0_n(1)+alpha1_n(1);
  gamma = 3*C2/(2*C1);
  alpha0(1) = alpha0_n(1)*gamma;
  alpha1(1) = alpha1_n(1)*gamma;
  alpha0(2) = alpha0_n(2)*gamma;
  alpha1(2) = alpha1_n(2)*gamma;
  bbeta0 = beta0_n*gamma;
  bbeta1 = beta1_n*gamma;
  bbeta2 = beta2_n*gamma;
 
 % Compute the betas, evaluating alphas and betas at 0, 1/2, 1.
  alpha_R_0 = alpha0(1);
  alpha_R_1 = alpha1(1);
  alpha_R_12 = (alpha0(1)+alpha1(1))/2;
  alpha_L_0 = alpha0(2);
  alpha_L_1 = alpha1(2);
  alpha_L_12 = (alpha0(2)+alpha1(2))/2;
  beta_0 = bbeta0;
  beta_1 = bbeta2;
  beta_12 = (bbeta0+bbeta2)/4+bbeta1/2;
 
 % Compute the matrix of the system considering the relationship between beta^(i,0), beta^(i,1) and beta
  M = [alpha_R_0 0 alpha_L_0 0; ...
       0 alpha_R_1 0 alpha_L_1; ...
       alpha_R_12/2 alpha_R_12/2 alpha_L_12/2 alpha_L_12/2];
 
  if (rank(M)==3)
 % Compute beta1_R, beta0_L, beta1_L in terms of beta0_R
    quant1 = (-alpha_L_12/2 + (alpha_L_0*alpha_R_12)/(2*alpha_R_0)) / ...
      (-(alpha_L_1*alpha_R_12)/(2*alpha_R_1) + alpha_L_12/2);
    quant2 = (beta_12-(beta_0*alpha_R_12)/(2*alpha_R_0) - (beta_1*alpha_R_12)/(2*alpha_R_1)) / ...
      (-(alpha_L_1*alpha_R_12)/(2*alpha_R_1) + alpha_L_12/2); 
 
 % beta1_R=a+b*beta0_R,  beta0_L=c+d*beta0_R,  beta1_L=e+f*beta0_R, where
    a = quant2; b = quant1;
    c = beta_0/alpha_R_0; d = -alpha_L_0/alpha_R_0;
    e = (beta_1-alpha_L_1*quant2)/alpha_R_1; f = -alpha_L_1*quant1/alpha_R_1;
      
 %We determine beta0_R by minimizing the sum of the norms of beta_R and beta_L
    C1 = ((b-1)^2)/3 + (b-1) + ((f-d)^2)/3 + (f-d)*d + d^2 + 1;
    C2 = 2*a*(b-1)/3 + a + 2*(e-c)*(f-d)/3 + (e-c)*d + (f-d)*c + 2*c*d;
    beta0(1) = -C2/(2*C1);
    beta1(1) = a + b*beta0(1);
    beta0(2) = c + d*beta0(1);
    beta1(2) = e + f*beta0(1);

  else
 % Compute beta0_L in terms of beta0_R and beta1_L in terms of beta1_R:
 % beta0_R=a+b*beta0_L,  beta1_R=c+d*beta1_L, where
    a = beta_0/alpha_R_0; b = -alpha_L_0/alpha_R_0;
    c = beta_1/alpha_R_1; d = -alpha_L_1/alpha_R_1;
 
 % Determine beta0_R and beta_1_R by minimizing the sum of the norms of beta_R and beta_L
    M2 = [2*(1+b^2) 1+b*d; 1+b*d 2*(1+d^2)];
    M2b = [-b*c-2*a*b; -a*d-2*c*d];
    sol = M2\M2b;
    beta0(1)= sol(1);
    beta1(1)= sol(2);
    beta0(2)= a + b*beta0(1);
    beta1(2)= c + d*beta1(1);
  end 
  
end


% Functions to deal with general orientation of the patches
function geo_reoriented = reorientation_patches (operations, geometry)
  nrb_patches = [geometry.nurbs];
  for ii = 1:numel(nrb_patches)
    if (operations(ii,1))
      nrb_patches(ii) = nrbreverse (nrb_patches(ii), 1);
    end
    if (operations(ii,2))
      nrb_patches(ii) = nrbreverse (nrb_patches(ii), 2);
    end
    if (operations(ii,3))
      nrb_patches(ii) = nrbtransp (nrb_patches(ii));
    end
  end

  geo_reoriented = mp_geo_load (nrb_patches);
%   [geo_reoriented, ~, local_interface] = mp_geo_load (nrb_patches);
%   if (numel (nrb_patches) == 2)
%     if (local_interface.side1 ~= 1 || local_interface.side2 ~= 3 || local_interface.ornt ~= 1)
%       error('The reorientation is wrong')
%     end
%   end
end

function knots = knot_vector_reorientation (knots, operations)
  if (operations(1))
    knots{1} = sort (1 - knots{1});
  end
  if (operations(2))
    knots{2} = sort (1-knots{2});
  end
  if (operations(3))
    knots([2 1]) = knots([1 2]);
  end
end

function degree = degree_reorientation (degree, transposition)
  if (transposition)
    degree = degree([2 1]);
  end
end

function indices = indices_reorientation (ndof_dir, operations)
  ndof = prod (ndof_dir);
  indices = reshape (1:ndof, ndof_dir);
  if (operations(1))
    indices = flipud (indices);
  end
  if (operations(2))
    indices = fliplr (indices);
  end
  if (operations(3))
    indices = indices.';
  end   
end
