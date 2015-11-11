function rhs = op_f_v_mp (space, msh, coeff, patch_list)

  if (nargin < 4)
    patch_list = 1:msh.npatch;
  end

  if (space.npatch ~= msh.npatch)
    error ('op_f_v_mp: the number of patches does not coincide')
  end
  
  rhs = zeros (space.ndof, 1);
  for iptc = patch_list
    rhs_loc = op_f_v_tp (space.sp_patch{iptc}, msh.msh_patch{iptc}, coeff);
    rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_loc;
  end

end