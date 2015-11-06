% TEST_SQUARE_G_NMNN: data function for Neumann boundary condition.

function g = test_square_g_nmnn (x, y, ind)
  switch ind
    case 1
      g = -exp (x) .* sin (y);
    case 2
      g = exp (x) .* sin (y);
    case 3
      g = -exp (x) .* cos (y);
    case 4
      g = exp (x) .* cos (y);
    otherwise
      error ('g_nmnn: unknown reference number');
  end
end

