% TEST_WAVEGUIDE_3D_H_DRCHLT: data function for Dirichlet boundary condition.

function h = test_waveguide_3d_h_drchlt (x, y, z, ind)

  switch (ind)
   case 3
    h = cat(1, reshape (0*x, [1, size(x)]), ...
            reshape (0*x+1, [1, size(x)]), ...
            reshape (0*x+2, [1, size(x)]));
    case 4
     h = cat(1, reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]));
    otherwise
     h = cat(1, reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]), ...
             reshape (0*x, [1, size(x)]));
  end

end

