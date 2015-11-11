function A = op_u_v_mp (spu, spv, msh, coeff, patch_list)

  if (nargin < 5)
    patch_list = 1:msh.npatch;
  end

  if ((spu.npatch ~= spv.npatch) || (spu.npatch ~= msh.npatch))
    error ('op_u_v_mp: the number of patches does not coincide')
  end
  
  ncounter = 0;
  for iptc = patch_list
    [rs, cs, vs] = op_u_v_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, coeff);
    rows(ncounter+(1:numel (rs))) = spv.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spu.gnum{iptc}(cs);
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spv.ndof, spu.ndof);
  clear rows cols vals rs cs vs

end