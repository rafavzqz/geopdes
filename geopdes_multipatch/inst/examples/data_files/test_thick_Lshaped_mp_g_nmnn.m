% TEST_THICL_LSHAPED_MP_G_NMNN: data function for Neumann boundary condition.

function g = test_thick_Lshaped_mp_g_nmnn(x, y, z, ind)
  switch (ind)
   case {1, 4}
     g = cos(z) .* exp(x) .* (sin(x.*y) + y .* cos(x.*y));
   case {2, 3}
     g = -x .* exp(x) .* cos(x.*y) .* cos(z);
   case {5}
     g = -cos(z) .* exp(x) .* (sin(x.*y) + y .* cos(x.*y));
   case {6}
     g = x .* exp(x) .* cos(x.*y) .* cos(z);
   case {7}
     g = exp(x) .* sin(x.*y) .* sin(z);
   case {8}
     g = - exp(x) .* sin(x.*y) .* sin(z);
   otherwise
    error('g_nmnn: error in the reference number for the boundary')
  end
end

