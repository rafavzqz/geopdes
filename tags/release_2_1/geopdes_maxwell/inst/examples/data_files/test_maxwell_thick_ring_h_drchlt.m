% TEST_MAXWELL_THICK_RING_H_DRCHLT: data function for Dirichlet boundary condition.

function h = test_maxwell_thick_ring_h_drchlt (x, y, z, ind)

  [theta, r] = cart2pol (x, y);
  h = zeros ([3, size(x)]);
  switch (ind)
    case 1
      h(1,:,:) = sin(theta) .* x .* y;
      h(2,:,:) = -cos(theta) .* x .* y;
      h(3,:,:) = -sin(theta) .* sin(y) + cos(theta) .* sin(x);
    case 2
      h(1,:,:) = -sin(theta) .* x .* y;
      h(2,:,:) = cos(theta) .* x .* y;
      h(3,:,:) = sin(theta) .* sin(y) - cos(theta) .* sin(x);
    case 3
      h(1,:,:) = x .* y;
      h(3,:,:) = -sin(y);
    case 4
      h(2,:,:) = -x .* y;
      h(3,:,:) = sin(x);
    case 5
      h(1,:,:) = -sin(x);
      h(2,:,:) = sin(y);
    case 6
      h(1,:,:) = sin(x);
      h(2,:,:) = -sin(y);
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end

