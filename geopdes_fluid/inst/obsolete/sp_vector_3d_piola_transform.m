% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_vector_3d_piola_transform (sp1, sp2, sp3, msh)

warning ('The class SP_VECTOR_3D_PIOLA_TRANSFORM is obsolete. Using SP_VECTOR_DIV_TRANSFORM instead. Read the help for the usage')

sp = sp_vector_div_transform ({sp1, sp2, sp3}, msh);
