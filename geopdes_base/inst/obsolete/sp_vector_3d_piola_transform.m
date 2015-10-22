% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_vector_3d_piola_transform (sp1, sp2, sp3, msh)

warning ('geopdes:obsolete', 'The class SP_VECTOR_3D_PIOLA_TRANSFORM is obsolete. Using SP_VECTOR with the transformation as the last argument. Read the help for the usage')

sp = sp_vector ({sp1, sp2, sp3}, msh, 'div-preserving');
