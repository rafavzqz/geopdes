% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_3d (interfaces, sp)

warning ('Function MP_INTERFACE_HCURL_3D is obsolete. Using MP_INTERFACE_HCURL instead')

[glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl (interfaces, sp);
