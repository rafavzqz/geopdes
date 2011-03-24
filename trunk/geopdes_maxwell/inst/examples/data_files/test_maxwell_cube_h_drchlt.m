% TEST_MAXWELL_CUBE_H_DRCHLT: data function for Dirichlet boundary condition.

function h = test_maxwell_cube_h_drchlt (x, y, z, ind)

  h = zeros ([3, size(x)]);
  switch (ind)
    case 1
      h(2,:,:) = -exp(z) .* cos(x);
      h(3,:,:) = exp(x) .* cos(y);
    case 2
      h(2,:,:) = exp(z) .* cos(x);
      h(3,:,:) = -exp(x) .* cos(y);
    case 3
      h(1,:,:) = exp(z) .* cos(x);
      h(3,:,:) = -sin(y) .* z;
    case 4
      h(1,:,:) = -exp(z) .* cos(x);
      h(3,:,:) = sin(y) .* z;
    case 5
      h(1,:,:) = -exp(x) .* cos(y);
      h(2,:,:) = sin(y) .* z;
    case 6
      h(1,:,:) = exp(x) .* cos(y);
      h(2,:,:) = -sin(y) .* z;
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end

