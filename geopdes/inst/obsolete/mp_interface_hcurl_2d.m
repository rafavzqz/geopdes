% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_2d (interfaces, sp)

warning ('geopdes:obsolete','Function MP_INTERFACE_HCURL_2D is obsolete. Using MP_INTERFACE_HCURL instead')

[glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl (interfaces, sp)
