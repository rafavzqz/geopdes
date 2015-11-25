% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = ...
              mp_solve_maxwell_eig_mixed1_2d (problem_data, method_data)

warning ('geopdes:obsolete','Function MP_SOLVE_MAXWELL_EIG_MIXED1_2D is obsolete. Using MP_SOLVE_MAXWELL_EIG_MIXED1 instead')

[geometry, msh, space, sp_mul, eigv, eigf, gnum, dofs_ornt, gnum_mul] = mp_solve_maxwell_eig_mixed1 (problem_data, method_data);

