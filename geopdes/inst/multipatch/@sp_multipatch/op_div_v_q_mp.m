function A = op_div_v_q_mp (spv, spq, msh, patch_list)

  if (nargin < 4)
    patch_list = 1:msh.npatch;
  end

  if ((spv.npatch ~= spq.npatch) || (spv.npatch ~= msh.npatch))
    error ('op_div_v_q_mp: the number of patches does not coincide')
  end
  
  ncounter = 0;
  for iptc = patch_list
    [rs, cs, vs] = op_div_v_q_tp (spv.sp_patch{iptc}, spq.sp_patch{iptc}, msh.msh_patch{iptc});
    rows(ncounter+(1:numel (rs))) = spq.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spv.gnum{iptc}(cs);

    if (~isempty (spv.dofs_ornt))
      vs = vs .* spv.dofs_ornt{iptc}(cs)';
    end
    
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spq.ndof, spv.ndof);
  clear rows cols vals rs cs vs

end