% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_vector_2d (sp1, sp2, msh)

warning ('geopdes:obsolete', 'The class SP_VECTOR_2D is obsolete. Using SP_VECTOR instead. Read the help for the usage')

sp = sp_vector ({sp1, sp2}, msh);
