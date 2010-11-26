% TEST_MAXWELL_SQUARE_G_NMNN: data function for Neumann boundary condition.

function g = test_maxwell_square_g_nmnn (x, y, ind)

  g = zeros ([2, size(x)]);
  switch (ind)
    case 1
      g(2,:,:) = -exp(x) .* cos(y) + cos(y);
    case 2
      g(2,:,:) = exp(x) .* cos(y) - cos(y);
    case 3
      g(1,:,:) = exp(x) .* cos(y) - cos(y);
    case 4
      g(1,:,:) = -exp(x) .* cos(y) + cos(y);
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end

