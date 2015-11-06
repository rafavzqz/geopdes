% TEST_THICK_RING_G_NMNN: data function for Neumann boundary condition.

function g = test_thick_ring_g_nmnn (x, y, z, ind)
  [theta, r] = cart2pol (x, y);
  switch ind
    case 1
      g = -cos (theta) .* exp (x) .* cos (z) .* (sin (x.*y) + y .* cos (x.*y)) -...
            sin (theta) .* x .* exp (x) .* cos (x.*y) .* cos (z);
    case 2
      g = cos (theta) .* exp (x) .* cos (z) .* (sin (x.*y) + y .* cos (x.*y)) +...
            sin (theta) .* x .* exp (x) .* cos (x.*y) .* cos (z);
    case 3
      g = -x .* exp (x) .* cos (x.*y) .* cos (z);
    case 4
      g = -exp (x) .* cos (z) .* (sin (x.*y) + y .* cos (x.*y));
    case 5
      g = exp (x) .* sin (x.*y) .* sin (z);
    case 6
      g = -exp (x) .* sin (x.*y) .* sin (z);
    otherwise
      error ('g_nmnn: unknown reference number');
  end
end
