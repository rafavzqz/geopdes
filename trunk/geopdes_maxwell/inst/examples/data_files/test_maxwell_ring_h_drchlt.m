% TEST_MAXWELL_RING_H_DRCHLT: data function for Dirichlet boundary condition.

function h = test_maxwell_ring_h_drchlt (x, y, ind)

  [theta, r] = cart2pol (x, y);
  h = zeros (size(x));
  switch (ind)
    case 1
      h = -sin(theta) .* sin(y) + cos(theta) .* sin(x);
    case 2
      h = sin(theta) .* sin(y) - cos(theta) .* sin(x);
    case 3
      h = -sin(y);
    case 4
      h = sin(x);
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end

