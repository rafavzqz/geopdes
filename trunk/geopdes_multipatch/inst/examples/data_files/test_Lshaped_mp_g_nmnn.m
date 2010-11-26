% TEST_LSHAPED_MP_G_NMNN: data function for Neumann boundary condition.

function g = test_Lshaped_mp_g_nmnn(x, y, ind, k)
  switch (ind)
   case {1, 6}
     g = exp(x) .* (sin (x.*y) + y .* cos(x.*y));
   case {2, 3}
     g = -x .* exp(x) .* cos(x.*y);
   case {4}
     g = -exp(x) .* (sin (x.*y) + y .* cos(x.*y));
   case {5}
     g = x .* exp(x) .* cos(x.*y);
   otherwise
    error('g_nmnn: error in the reference number for the boundary')
  end
end

