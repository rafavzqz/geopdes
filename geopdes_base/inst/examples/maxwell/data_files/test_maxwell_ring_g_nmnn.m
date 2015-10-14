% TEST_MAXWELL_RING_G_NMNN: data function for Neumann boundary condition.

function g = test_maxwell_ring_g_nmnn (x, y, ind)

  [theta, r] = cart2pol (x, y);
  g = zeros ([2, size(x)]);
  switch (ind)
    case 1
      g(1,:,:) = sin(theta) .* (cos(x) - cos(y));
      g(2,:,:) = cos(theta) .* (cos(y) - cos(x));
    case 2
      g(1,:,:) = sin(theta) .* (cos(y) - cos(x));
      g(2,:,:) = cos(theta) .* (cos(x) - cos(y));
    case 3
      g(1,:,:) = cos(x) - cos(y);
    case 4
      g(2,:,:) = cos(y) - cos(x);
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end

