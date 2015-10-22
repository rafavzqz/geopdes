% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_bspline_nforms (knots, degree, msh)

  warning ('geopdes:obsolete', 'The function SP_BSPLINE_NFORMS is obsolete. Using SP_BSPLINE with the transformation as last argument')
  sp = sp_bspline (knots, degree, msh, 'integral-preserving');

end
