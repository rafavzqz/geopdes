% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space, u] = ...
              solve_maxwell_src_3d (problem_data, method_data)

warning ('geopdes:obsolete', ['Function SOLVE_MAXWELL_SRC_3D is obsolete. Using SOLVE_MAXWELL_SRC instead. \n', ...
    'Be aware that the way to impose non-homegeneous Dirichlet BC has also changed.'])

[geometry, msh, space, u] = solve_maxwell_src (problem_data, method_data);

