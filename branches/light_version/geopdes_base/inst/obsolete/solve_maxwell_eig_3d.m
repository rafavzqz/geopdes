% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, eigv, eigf] = ...
              solve_maxwell_eig_3d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_MAXWELL_EIG_3D is obsolete. Using SOLVE_MAXWELL_EIG instead')

[geometry, msh, space, eigv, eigf] = solve_maxwell_eig (problem_data, method_data);

