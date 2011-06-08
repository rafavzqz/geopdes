% RING_POLAR_MAP_DER: Jacobian of the parameterization given in ring_polar_map.

function jac = ring_polar_map_der (pts)
  u = pts(1,:); v = pts(2,:);
  jac = zeros (2, 2, numel (u));
  jac(1,1,:) =  cos (pi.*v/2);
  jac(2,1,:) =  sin (pi.*v/2);
  jac(1,2,:) = -pi*(u+1) .* sin (pi.*v/2)/2;
  jac(2,2,:) =  pi*(u+1) .* cos (pi.*v/2)/2;
end

