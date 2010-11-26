function uex = test_plane_strain_ring_uex (x, y, E, nu, P)

% Inner and outer radius of our ring
  R_i = 1; R_o = 2;

  [theta, rad] = cart2pol (x, y);
  theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
  
  u_r = P*R_i^2/(E*(R_o^2 - R_i^2)) * ((1-nu)*rad + (1+nu)*R_o^2./rad);
  uex = zeros (2, size (x, 1), size (x, 2));
  uex(1, :, :) = u_r .* cos (theta);
  uex(2, :, :) = u_r .* sin (theta);

end
