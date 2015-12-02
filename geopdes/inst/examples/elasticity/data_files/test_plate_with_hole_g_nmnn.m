function g = test_plate_with_hole_g_nmnn (x, y, ind)

% Inner radius and applied stress
  Ri = 1; Tx = 10;

  switch ind
    case 4
      g_polar = zeros(2, size(x,1), size(x,2));
      g = zeros(2, size(x,1), size(x,2));
      [theta, rad] = cart2pol (x, y);

      S_rr = Tx/2 * ((1 - (Ri./rad).^2) + cos(2*theta) .* (1 - 4*(Ri./rad).^2 + 3*(Ri./rad).^4));
      S_tt = Tx/2 * ((1 + (Ri./rad).^2) - cos(2*theta) .* (1 + 3*(Ri./rad).^4));
      S_rt = -Tx/2 * sin (2*theta) .* (1 + 2*(Ri./rad).^2 - 3*(Ri./rad).^4);

      g_polar(1,:,:) = (-S_rr .* cos(theta) + S_rt .* sin(theta)) .* (theta > 3*pi/4) +  ...
                        (S_rr .* sin(theta) + S_rt .* cos(theta)) .* (theta < 3*pi/4);
      g_polar(2,:,:) = (-S_rt .* cos(theta) + S_tt .* sin(theta)) .* (theta > 3*pi/4) + ...
                        (S_rt .* sin(theta) + S_tt .* cos(theta)) .* (theta < 3*pi/4);
      g(1,:,:) = reshape(g_polar(1,:,:), size(x)) .* cos(theta) - reshape(g_polar(2,:,:), size(x)) .* sin(theta);
      g(2,:,:) = reshape(g_polar(1,:,:), size(x)) .* sin(theta) + reshape(g_polar(2,:,:), size(x)) .* cos(theta);
    otherwise
      error ('g_nmnn: unknown reference number');
  end
end

