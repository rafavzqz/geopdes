function varargout = op_su_ev_tp (space1, space2, msh, lambda, mu)

  A = spalloc (space2.ndof, space1.ndof, 5*space1.ndof);

  ndim = numel (msh.qn);

  for iel = 1:msh.nelu
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, ...
                               'gradient', true, 'divergence', true);

    for idim = 1:ndim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A = A + op_su_ev (sp1_col, sp2_col, msh_col, lambda (x{:}), mu (x{:}));
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
