% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_vector_2d_piola_transform (sp1, sp2, msh)

warning ('geopdes:obsolete', 'The class SP_VECTOR_2D_PIOLA_TRANSFORM is obsolete. Using SP_VECTOR with the transformation as last argument')

sp = sp_vector ({sp1, sp2}, msh, 'div-preserving');
