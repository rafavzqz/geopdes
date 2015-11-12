function A = op_v_gradp_mp (spv, spp, msh, coeff, patch_list)

  if (nargin < 5)
    patch_list = 1:msh.npatch;
  end

  if ((spv.npatch ~= spp.npatch) || (spv.npatch ~= msh.npatch))
    error ('op_v_grad_mp: the number of patches does not coincide')
  end
  
  ncounter = 0;
  for iptc = patch_list
    [rs, cs, vs] = op_v_gradp_tp (spv.sp_patch{iptc}, spp.sp_patch{iptc}, msh.msh_patch{iptc}, coeff);
    rows(ncounter+(1:numel (rs))) = spp.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spv.gnum{iptc}(cs);

    if (~isempty (spv.dofs_ornt))
      vs = vs .* spv.dofs_ornt{iptc}(cs)';
    end
    
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spp.ndof, spv.ndof);
  clear rows cols vals rs cs vs

end