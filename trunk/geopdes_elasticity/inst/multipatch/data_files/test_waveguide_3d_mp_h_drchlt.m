% TEST_WAVEGUIDE_3D_MP_H_DRCHLT: data function for Dirichlet boundary condition.

function h = test_waveguide_3d_mp_h_drchlt (x, y, z, ind)
  switch (ind)
   case {1, 11, 12, 13, 14}
    h = cat(1, reshape (0*x, [1, size(x)]), ...
            reshape (0*x+.5, [1, size(x)]), ...
            reshape (0*x+.7, [1, size(x)]));
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

