% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, sp, u, gnum] = ...
                    mp_solve_plane_strain_2d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function MP_SOLVE_PLANE_STRAIN_2D is obsolete. Using MP_SOLVE_LINEAR_ELASTICITY instead')

[geometry, msh, sp, u, gnum] = mp_solve_linear_elasticity (problem_data, method_data);

