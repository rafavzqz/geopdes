% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [spv, spp, PI] = sp_bspline_fluid_2d (element_name, ...
                   knots, nsub_p, degree_p, regularity_p, msh)

warning ('geopdes:obsolete', 'Function SP_BSPLINE_FLUID_2D is obsolete. Using SP_BSPLINE_FLUID instead')

[spv, spp] = sp_bspline_fluid (element_name, knots, nsub_p, degree_p, regularity_p, msh);

if (nargout == 3)
  switch (lower (element_name))
    case {'th', 'sg', 'ndl'}
      PI = speye (spp.ndof);
    case {'rt'}
     warning ('The use of the T-spline space for the pressure with RT-splines is deprecated. Try imposing Dirichlet boundary conditions in a weak form (see solve_stokes)')
      PI = b2nst__ (spp, spp.knots, degree_p, msh);
  end
end
