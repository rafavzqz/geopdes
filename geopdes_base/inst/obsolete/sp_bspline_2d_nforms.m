% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_bspline_2d_nforms (knots, degree, msh)

warning ('The class SP_BSPLINE_2D_nforms is obsolete. Using SP_BSPLINE_NFORMS instead')

sp = sp_bspline_nforms (knots, degree, msh);
