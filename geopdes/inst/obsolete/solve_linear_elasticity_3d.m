% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, sp, u] = solve_linear_elasticity_3d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_LINEAR_ELASTICITY_3D is obsolete. Using SOLVE_LINEAR_ELASTICITY instead')

[geometry, msh, sp, u] = solve_linear_elasticity (problem_data, method_data);

