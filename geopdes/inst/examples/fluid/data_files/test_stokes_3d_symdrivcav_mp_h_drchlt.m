% TEST_STOKES_3D_SYMDRIVCAV_MP_H_DRCHLT: data function for Dirichlet boundary condition.

function h  = test_stokes_3d_symdrivcav_mp_h_drchlt (x, y, z, iside) 
switch (iside)
 case {1, 2, 5, 6}
  h = zeros ([3, size(x)]);
 case 3
  h1 = -ones ([1, size(x)]);
  h2 = zeros ([1, size(x)]);
  h3 = zeros ([1, size(x)]);
  h  = cat (1, h1, h2, h3);
 case 4
  h1 = ones ([1, size(x)]);
  h2 = zeros ([1, size(x)]);
  h3 = zeros ([1, size(x)]);
  h  = cat (1, h1, h2, h3);
end
end
