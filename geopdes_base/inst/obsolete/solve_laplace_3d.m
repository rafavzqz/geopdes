% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, u] = ...
              solve_laplace_3d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_LAPLACE_3D is obsolete. Using SOLVE_LAPLACE instead')

[geometry, msh, space, u] = solve_laplace (problem_data, method_data);

