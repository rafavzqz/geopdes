% TEST_RING_MIXED_BC_G_NMNN: data function for Neumann boundary condition.

function g = test_ring_mixed_bc_g_nmnn (x, y, ind)

  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      g = -cos (theta) .* exp (x) .* (sin (x.*y) + y.*cos (x.*y)) - ...
          sin (theta) .* exp(x) .* x .* cos (x.*y);
    case 2
      g = cos (theta) .* exp (x) .* (sin (x.*y) + y.*cos (x.*y)) + ...
          sin (theta) .* exp(x) .* x .* cos (x.*y);
    case 3
      g = -x .* exp(x) .* cos (x.*y);
    case 4
      g = -exp(x) .* (sin (x.*y) + y .* cos (x.*y));
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end

