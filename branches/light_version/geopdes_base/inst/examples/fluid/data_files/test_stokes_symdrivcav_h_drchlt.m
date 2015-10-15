% TEST_STOKES_SYMDRIVCAV_H_DRCHLT: data function for Dirichlet boundary condition.

function h  = test_stokes_symdrivcav_h_drchlt (x, y, iside) 
switch (iside)
 case {1, 2}
  h = zeros ([2, size(x)]);
 case 3
  h1 = -ones ([1, size(x)]);
  h2 = zeros ([1, size(x)]);
  h = cat (1, h1, h2);
 case 4
  h1 = ones ([1, size(x)]);
  h2 = zeros ([1, size(x)]);
  h = cat (1, h1, h2);
end
end
