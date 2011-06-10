function varargout = op_su_ev_tp (space1, space2, msh, lambda, mu)

  A = spalloc (space2.ndof, space1.ndof, 5*space1.ndof);

  for iel = 1:msh.nelu
    [sp1, element_list] = sp_evaluate_col (space1, msh, iel, 'value', false, ...
                                        'gradient', true, 'divergence', true);
    sp2 = sp_evaluate_col (space2, msh, iel, 'value', false, ... 
                                        'gradient', true, 'divergence', true);

    A = A + op_su_ev (sp1, sp2, msh, lambda, mu, element_list);
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
