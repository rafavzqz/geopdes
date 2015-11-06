% TEST_MAXWELL_CUBE_G_NMNN: data function for Neumann boundary condition.

function g = test_maxwell_cube_g_nmnn (x, y, z, ind)

  g = zeros ([3, size(x)]);
  switch (ind)
    case 1
      g(2,:,:) = -(exp(x) .* cos(y) - z .* cos(y));
      g(3,:,:) = sin(y) + exp(z) .* sin(x);
    case 2
      g(2,:,:) = exp(x) .* cos(y) - z .* cos(y);
      g(3,:,:) = -(sin(y) + exp(z) .* sin(x));
    case 3
      g(1,:,:) = exp(x) .* cos(y) - z .* cos(y);
    case 4
      g(1,:,:) = -(exp(x) .* cos(y) - z .* cos(y));
    case 5
      g(1,:,:) = -(sin(y) + exp(z) .* sin(x));
    case 6
      g(1,:,:) = sin(y) + exp(z) .* sin(x);
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end

