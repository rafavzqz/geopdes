% TEST_NAVIERSTOKES_SQUARE_GRADUEX: gradient of the exact solution.

function gu = test_navierstokes_square_graduex (x, y)

  uxx = @(x, y) 0*x;
  uxy = @(x, y) -pi*sin(y*pi);
  uyx = @(x, y) 2*x - x.^0;
  uyy = @(x, y) 0*x;

  gu = zeros (2, 2, size(x,1), size(x,2));
  gu(1, 1, :, :) = reshape (uxx (x,y), [1, 1, size(x)]);
  gu(2, 1, :, :) = reshape (uyx (x,y), [1, 1, size(x)]);
  gu(1, 2, :, :) = reshape (uxy (x,y), [1, 1, size(x)]);
  gu(2, 2, :, :) = reshape (uyy (x,y), [1, 1, size(x)]);
end
