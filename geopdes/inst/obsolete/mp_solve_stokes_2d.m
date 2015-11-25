% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [geometry, msh, spv, vel, gnum, spp, press, gnump] = ...
              mp_solve_stokes_2d (problem_data, method_data)

warning ('geopdes:obsolete','Function MP_SOLVE_STOKES_2D is obsolete. Using MP_SOLVE_STOKES instead')

function [geometry, msh, spv, vel, gnum, spp, press, gnump] = ...
              mp_solve_stokes (problem_data, method_data)

