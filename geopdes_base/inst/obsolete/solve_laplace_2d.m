% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, u] = ...
              solve_laplace_2d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_LAPLACE_2D is obsolete. Using SOLVE_LAPLACE instead')

[geometry, msh, space, u] = solve_laplace (problem_data, method_data);

