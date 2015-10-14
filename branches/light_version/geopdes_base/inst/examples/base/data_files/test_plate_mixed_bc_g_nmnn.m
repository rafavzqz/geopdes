% TEST_PLATE_MIXED_BC_G_NMNN: data function for Neumann boundary condition.

function g = test_plate_mixed_bc_g_nmnn (x, y, ind)
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      g = -exp(x) .* cos (y);
    case 2
      g = exp(x) .* sin (y);
    case 3
      g = -cos (theta) .* exp (x) .* sin (y) - ...
          sin (theta) .* exp(x) .* cos (y);
    case 4
      g = (theta > 3*pi/4) .* (-exp(x) .* sin (y))  + ...
          (theta < 3*pi/4) .* (exp(x) .* cos (y));
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end

