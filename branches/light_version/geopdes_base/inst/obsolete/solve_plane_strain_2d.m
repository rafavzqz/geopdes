% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, sp, u] = solve_plane_strain_2d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_PLANE_STRAIN_2D is obsolete. Using SOLVE_LINEAR_ELASTICITY instead')

[geometry, msh, sp, u] = solve_linear_elasticity (problem_data, method_data);

