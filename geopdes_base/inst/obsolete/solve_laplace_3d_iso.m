% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, u] = ...
              solve_laplace_3d_iso (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_LAPLACE_3D_ISO is obsolete. Using SOLVE_LAPLACE_ISO instead')

[geometry, msh, space, u] = solve_laplace_iso (problem_data, method_data);
