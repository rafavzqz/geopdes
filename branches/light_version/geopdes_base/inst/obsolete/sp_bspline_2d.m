% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_bspline_2d (knots, degree, msh)

warning ('geopdes:obsolete', 'The class SP_BSPLINE_2D is obsolete. Using SP_BSPLINE instead')

sp = sp_bspline (knots, degree, msh);
