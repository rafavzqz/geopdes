% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_bspline_2d_nforms (knots, degree, msh)

warning ('geopdes:obsolete', 'The class SP_BSPLINE_2D_nforms is obsolete. Using SP_BSPLINE with the transformation as last argument')

sp = sp_bspline_nforms (knots, degree, msh, 'integral-preserving');
