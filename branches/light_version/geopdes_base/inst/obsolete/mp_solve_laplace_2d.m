% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, sp, u, gnum] = ...
              mp_solve_laplace_2d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function MP_SOLVE_LAPLACE_2D is obsolete. Using MP_SOLVE_LAPLACE instead')

[geometry, msh, sp, u, gnum] = mp_solve_laplace (problem_data, method_data);
