% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, eigv, eigf, gnum, dofs_ornt] = ...
              mp_solve_maxwell_eig_3d (problem_data, method_data)

warning ('Function MP_SOLVE_MAXWELL_EIG_3D is obsolete. Using SOLVE_MAXWELL_EIG instead')

[geometry, msh, space, eigv, eigf, gnum, dofs_ornt] = mp_solve_maxwell_eig (problem_data, method_data);

