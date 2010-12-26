% TEST_SQUARE_HEAT_EXP_G_NMNN: data function for Neumann boundary condition.

function g = test_square_heat_exp_g_nmnn (x, y, t, ind)
  switch ind
    case 1
      g = -2 * x .* exp(y) * exp(t);
    case 2
      g = 2 * x .* exp(y) * exp(t);
    case 3
      g = -x .* x .* exp(y) * exp(t);
    case 4
      g = x .* x .* exp(y) * exp(t);
    otherwise
      error ('g_nmnn: unknown reference number');
  end
end

