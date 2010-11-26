% TEST_STOKES_SQUARE_BC_GRADUEX: gradient of the exact solution.

function gu = test_stokes_square_bc_graduex (x, y)
  uxx = @(x, y) (cos(x) .* cos(y));
  uxy = @(x, y) (-sin(x) .* sin(y));
  uyx = @(x, y) (sin(x) .* sin(y));
  uyy = @(x, y) (-cos(x) .* cos(y));
  
  gu = zeros (2, 2, size(x,1), size(x,2));
  gu(1, 1, :, :) = reshape (uxx (x,y), [1, 1, size(x)]);
  gu(2, 1, :, :) = reshape (uyx (x,y), [1, 1, size(x)]);
  gu(1, 2, :, :) = reshape (uxy (x,y), [1, 1, size(x)]);
  gu(2, 2, :, :) = reshape (uyy (x,y), [1, 1, size(x)]);
end
