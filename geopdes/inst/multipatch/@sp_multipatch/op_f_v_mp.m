function rhs = op_f_v_mp (space, msh, coeff)

  if (space.npatch ~= msh.npatch)
    error ('op_gradu_gradv_mp: the number of patches does not coincide')
  end
  
  rhs = zeros (space.ndof, 1);
  for iptc = 1:msh.npatch
    rhs_loc = op_f_v_tp (space.sp_patch{iptc}, msh.msh_patch{iptc}, coeff);
    rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_loc;
  end

end