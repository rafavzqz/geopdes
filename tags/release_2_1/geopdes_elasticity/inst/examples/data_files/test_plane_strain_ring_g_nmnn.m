function g = test_plane_strain_ring_g_nmnn (x, y, P, nu, ind)

% Inner and outer radius of our ring
  R_i = 1; R_o = 2;

  switch ind
    case 1
      g = zeros(2, size(x,1), size(x,2));
      [theta, rad] = cart2pol (x, y);
      theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;

      C = (P * R_i^2) / (R_o^2 - R_i^2);
      C = C*((1-nu)/((1+nu)*(1-2*nu)) - R_o^2/R_i^2);

      g(1, :, :) = -C * cos (theta);
      g(2, :, :) = -C * sin (theta);

    case 2
      g = zeros(2, size(x,1), size(x,2));
      [theta, rad] = cart2pol (x, y);
      theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;

      C = (P * R_i^2) / (R_o^2 - R_i^2);
      C = C*((1-nu)/((1+nu)*(1-2*nu)) - 1);

      g(1, :, :) = C * cos (theta);
      g(2, :, :) = C * sin (theta);

    otherwise
      error ('g_nmnn: unknown reference number');
  end
end

