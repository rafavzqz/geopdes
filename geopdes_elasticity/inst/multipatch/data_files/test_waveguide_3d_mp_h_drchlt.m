% TEST_WAVEGUIDE_3D_MP_H_DRCHLT: data function for Dirichlet boundary condition.

function h = test_waveguide_3d_mp_h_drchlt (x, y, z, ind)

  a = cos (pi/8);
  b = sin (pi/8);
  switch (ind)
   case {1, 11, 12, 13, 14}
    h = cat(1, reshape (0*x, [1, size(x)]), ...
            reshape ((a-1)*y-b*z, [1, size(x)]), ...
            reshape (b*y+(a-1)*z, [1, size(x)]));
    case {2, 7, 8, 9, 10}
     h = cat(1, reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]));
    otherwise
     h = cat(1, reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]));
  end

end

