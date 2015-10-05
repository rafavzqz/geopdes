% TEST_MAXWELL_SQUARE_H_DRCHLT: data function for Dirichlet boundary condition.

function h = test_maxwell_square_h_drchlt (x, y, ind)

  h = zeros (size (x));
  switch (ind)
    case 1
      h = exp(x) .* cos(y);
    case 2
      h = -exp(x) .* cos(y);
    case 3
      h = -sin(y);
    case 4
      h = sin(y);
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end

