% TEST_STOKES_BIFURCATION_MP_H_DRCHLT: data function for Dirichlet boundary condition.

function h  = test_stokes_bifurcation_mp_h_drchlt (x, y, iside) 
  switch (iside)
    case 2

      h = zeros ([2, size(x)]);

    case 1

      h2 = zeros (size(x));
      h1 = 1 - (y/.1).^2;
      h =cat(1, reshape (h1, [1, size(x)]), reshape (h2, [1, size(x)]));

    case 3

      h1 = -100*y.^2 + 160 * abs (y) - 63;
      h2 = zeros (size(x));
      h =cat(1, reshape (h1, [1, size(x)])/2, reshape (h2, [1, size(x)]));

  end
end
