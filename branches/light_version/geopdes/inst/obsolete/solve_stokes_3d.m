% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, space_v, vel, space_p, press] = ...
                          solve_stokes_3d (problem_data, method_data)

warning ('geopdes:obsolete', 'Function SOLVE_STOKES_3D is obsolete. Using SOLVE_STOKES instead')

[geometry, msh, space_v, vel, space_p, press] = solve_stokes (problem_data, method_data);

