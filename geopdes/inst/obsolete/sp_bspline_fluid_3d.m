% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [spv, spp] = sp_bspline_fluid_3d (element_name, ...
                   knots, nsub_p, degree_p, regularity_p, msh)

warning ('geopdes:obsolete', 'Function SP_BSPLINE_FLUID_3D is obsolete. Using SP_BSPLINE_FLUID instead')

[spv, spp] = sp_bspline_fluid (element_name, knots, nsub_p, degree_p, regularity_p, msh);
