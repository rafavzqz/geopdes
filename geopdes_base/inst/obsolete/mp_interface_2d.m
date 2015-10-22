% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [glob_num, glob_ndof] = mp_interface_2d (interfaces, sp)

warning ('geopdes:obsolete', 'Function MP_INTERFACE_2D is obsolete. Using MP_INTERFACE instead')

[glob_num, glob_ndof] = mp_interface (interfaces, sp);
