% TEST_CUBE_G_NMNN: data function for Neumann boundary condition.

function g = test_cube_g_nmnn (x, y, z, ind)
  switch ind
    case 1
      g = -exp (x + z) .* sin (y);
    case 2
      g = exp (x + z) .* sin (y);
    case 3
      g = -exp (x + z) .* cos (y);
    case 4
      g = exp (x + z) .* cos (y);
    case 5
      g = -exp (x + z) .* sin (y);
    case 6
      g = exp (x + z) .* sin (y);
    otherwise
      error ('test_cube_g_nmnn: unknown reference number');
  end
end
