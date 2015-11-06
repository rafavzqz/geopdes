% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [glob_num, glob_ndof] = mp_interface_vector_2d (interfaces, sp)

warning ('geopdes:obsolete', 'Function MP_INTERFACE_VECTOR_2D is obsolete. Using MP_INTERFACE_VECTOR instead')

[glob_num, glob_ndof] = mp_interface_vector (interfaces, sp);
