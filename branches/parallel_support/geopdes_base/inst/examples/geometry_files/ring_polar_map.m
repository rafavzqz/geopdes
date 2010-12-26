% RING_POLAR_MAP: parameterization of a ring using polar coordinates.

function F = ring_polar_map (pts)
  u = pts(1,:); v = pts(2,:);
  F(1,:) = (u+1) .* cos (pi.*v/2);
  F(2,:) = (u+1) .* sin (pi.*v/2);
end

