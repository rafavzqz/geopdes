% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, u] = ...
              solve_maxwell_src_2d (problem_data, method_data)

warning ('Function SOLVE_MAXWELL_SRC_2D is obsolete. Using SOLVE_MAXWELL_SRC instead')

[geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);

