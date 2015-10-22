% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, lambda, u] = ...
              solve_laplace_eig_2d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_LAPLACE_EIG_2D is obsolete. Using SOLVE_LAPLACE_EIG instead')

[geometry, msh, space, lambda, u] = solve_laplace_eig (problem_data, method_data);
