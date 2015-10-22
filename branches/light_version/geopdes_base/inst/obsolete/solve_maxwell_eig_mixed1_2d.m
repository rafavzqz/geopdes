% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, sp_mul, eigv, eigf] = ...
              solve_maxwell_eig_mixed1_2d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_MAXWELL_EIG_MIXED1_2D is obsolete. Using SOLVE_MAXWELL_EIG_MIXED1 instead')

[geometry, msh, space, sp_mul, eigv, eigf] = solve_maxwell_eig_mixed1 (problem_data, method_data);

