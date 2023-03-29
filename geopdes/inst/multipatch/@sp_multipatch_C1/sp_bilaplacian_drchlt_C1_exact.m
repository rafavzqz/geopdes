function [u_drchlt, drchlt_dofs, kernel_info] = sp_bilaplacian_drchlt_C1_exact_rafa_mp (space, msh, refs, uex, gradex)

% refs should be the whole boundary, for now
M = spalloc (space.ndof, space.ndof, space.ndof);
rhs = zeros (space.ndof, 1);

M2 = spalloc (space.ndof, space.ndof, space.ndof);
rhs2 = zeros (space.ndof, 1);

drchlt_dofs = [];
drchlt_dofs2 = [];

boundaries = msh.boundaries;
for iref = refs
%   href = @(varargin) h(varargin{:}, iref);
  for bnd_side = 1:boundaries(iref).nsides
    iptc = boundaries(iref).patches(bnd_side);
    iside = boundaries(iref).faces(bnd_side);

    msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside);
    msh_side_from_interior = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside);
    sp_bnd = space.sp_patch{iptc}.constructor (msh_side_from_interior);
    sp_bnd_struct = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);

    [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
    [~,icol] = find (Cpatch(sp_bnd.boundary(iside).dofs,:));
    [~,jcol] = find (Cpatch(sp_bnd.boundary(iside).adjacent_dofs,:));
    
    drchlt_dofs = union (drchlt_dofs, Cpatch_cols(icol));
    drchlt_dofs2 = union (drchlt_dofs2, Cpatch_cols(jcol));
    
    for idim = 1:msh.rdim
      x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    coeff_at_qnodes = ones (size(x{1}));
    dudn_at_qnodes = reshape (sum (gradex(x{:}) .* msh_side.normal, 1), msh_side.nqn, msh_side.nel) .* msh_side.charlen;

    M(Cpatch_cols,Cpatch_cols) = M(Cpatch_cols, Cpatch_cols) + ...
      Cpatch.' * op_u_v (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes) * Cpatch;
    rhs(Cpatch_cols) = rhs(Cpatch_cols) + Cpatch.' * op_f_v (sp_bnd_struct, msh_side, uex(x{:}));
    
    M2(Cpatch_cols,Cpatch_cols) = M2(Cpatch_cols,Cpatch_cols) + ...
      Cpatch.' * op_gradu_n_gradv_n (sp_bnd_struct, sp_bnd_struct, msh_side, coeff_at_qnodes.*msh_side.charlen) * Cpatch;
    rhs2(Cpatch_cols) = rhs2(Cpatch_cols) + Cpatch.' * op_gradv_n_f (sp_bnd_struct, msh_side, dudn_at_qnodes); % I am missing the other part of the vector. It is in M2 :-)    
  end
end

M_bdry = M + M2;
dofs_to_remove = [];
vertices_numbers = [];
row_indices = [];
count_vert = 0;

% Check the kernel of vertex functions on Dirichlet boundary vertices
% Pick up the basis function with the max. abs. coeff in the kernel, 
%  remove it from drchlt_dofs, and add the function in the kernel into the
%  internal part (it goes in the output)
B_change_local = [];
n_boundaries = numel(msh.boundaries); % number of boundary edges
global_refs = numel(space.interfaces) - n_boundaries + refs; % global numbering of Dirichlet boundary edges

for iv = 1 : numel(space.vertices)
  if (space.vertices(iv).boundary_vertex && space.vertices(iv).valence_p>1)
  % Loop just over Dirichlet boundary vertices
    if (~isempty(intersect(global_refs, space.vertices(iv).edges)))
      patches = space.vertices(iv).patches([1 end]);

      operations = space.vertices(iv).patch_reorientation([1 end], :);
      indices_loc_R = indices_reorientation(space.sp_patch{patches(1)}.ndof_dir, operations(1, :));
      indices_loc_L = indices_reorientation(space.sp_patch{patches(2)}.ndof_dir, operations(2, :));

      indices_loc_R = indices_loc_R(:);
      indices_loc_L = indices_loc_L(:);

%       Cpatch_ind_R = indices_loc_R([2 3 space.sp_patch{patches(1)}.ndof_dir(1)+2]);
%       Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(2)}.ndof_dir(1)+1 space.sp_patch{patches(2)}.ndof_dir(1)+2 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
      Cpatch_ind_R = indices_loc_R([1 2 3 space.sp_patch{patches(1)}.ndof_dir(1)+[1 2]]);
      if (space.vertices(iv).valence_p == 2)
        Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
      else
        Cpatch_ind_L = indices_loc_L([space.sp_patch{patches(2)}.ndof_dir(1)+[1 2] 2*space.sp_patch{patches(2)}.ndof_dir(1)+1]);
      end

% TODO: it is probably enough to use CC_vertices (and CC_edges)
      [Cpatch1, Cpatch_cols1] = sp_compute_Cpatch (space, patches(1));
      [Cpatch2, Cpatch_cols2] = sp_compute_Cpatch (space, patches(2));

      [~,~,inds1] = intersect (space.dofs_on_vertex{iv}, Cpatch_cols1);
      [~,~,inds2] = intersect (space.dofs_on_vertex{iv}, Cpatch_cols2);

      M_ker = [Cpatch1(Cpatch_ind_R, inds1); ...
               Cpatch2(Cpatch_ind_L, inds2)];
%       M_ker = M_bdry(space.dofs_on_vertex{iv}, space.dofs_on_vertex{iv});
      ker = null(full(M_ker));
      if (~isempty(ker))
        count_vert = count_vert + 1;
        [~, ind] = max(abs(ker));

        row_inds = (count_vert-1)*6 + (1:6);
        B_change_local = blkdiag (B_change_local, ker(:));

        dofs_on_vertex = space.dofs_on_vertex{iv};
        vertices_numbers(count_vert) = iv;
        dofs_to_remove(count_vert) = dofs_on_vertex(ind);
        row_indices(row_inds) = dofs_on_vertex;
      end
    end
  end
end

kernel_info = struct ('vertices_numbers', vertices_numbers, 'all_vertex_dofs', row_indices, 'quasi_interior_dofs', dofs_to_remove, 'B_change_local', sparse(B_change_local));

drchlt_dofs = union (drchlt_dofs, drchlt_dofs2);
drchlt_dofs = setdiff(drchlt_dofs, dofs_to_remove);

u_drchlt = M_bdry(drchlt_dofs,drchlt_dofs) \ (rhs(drchlt_dofs) + rhs2(drchlt_dofs));

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
