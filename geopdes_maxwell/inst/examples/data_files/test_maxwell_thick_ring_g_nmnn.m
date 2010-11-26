% TEST_MAXWELL_THICK_RING_G_NMNN: data function for Neumann boundary condition.

function g = test_maxwell_thick_ring_g_nmnn (x, y, z, ind)

  [theta, r] = cart2pol (x, y);
  g = zeros ([3, size(x)]);
  switch (ind)
    case 1
      g(1,:,:) = sin(theta) .* (cos(x) - cos(y));
      g(2,:,:) = cos(theta) .* (cos(y) - cos(x));
      g(3,:,:) = -x .* sin(theta) - y .* cos(theta);
    case 2
      g(1,:,:) = sin(theta) .* (cos(y) - cos(x));
      g(2,:,:) = cos(theta) .* (cos(x) - cos(y));
      g(3,:,:) = x .* sin(theta) + y .* cos(theta);
    case 3
      g(1,:,:) = cos(x) - cos(y);
    case 4
      g(2,:,:) = cos(y) - cos(x);
    case 5
      g(1,:,:) = y;
      g(2,:,:) = x;
    case 6
      g(1,:,:) = -y;
      g(2,:,:) = -x;
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end

