% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_bspline_3d (knots, degree, msh)

warning ('The class SP_BSPLINE_3D is obsolete. Using SP_BSPLINE instead')

sp = sp_bspline (knots, degree, msh);
