% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_vector_3d_curl_transform (sp1, sp2, sp3, msh)

warning ('geopdes:obsolete','The class SP_VECTOR_3D_CURL_TRANSFORM is obsolete. Using SP_VECTOR_CURL_TRANSFORM instead. Read the help for the usage')

sp = sp_vector_curl_transform ({sp1, sp2, sp3}, msh);
