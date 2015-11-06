% RING_POLAR_MAP_DER: Jacobian of the parameterization given in ring_polar_map.

function varargout = ring_polar_map_der (pts)

  if (iscell (pts))
    u = reshape (repmat (pts{1}(:), 1, numel (pts{2})), 1, []);
    v = reshape (repmat (pts{2}(:)', numel (pts{1}), 1), 1, []);
  else
    u = pts(1,:); v = pts(2,:);
  end
  
  jac = zeros (2, 2, numel (u));
  jac(1,1,:) =  cos (pi.*v/2);
  jac(2,1,:) =  sin (pi.*v/2);
  jac(1,2,:) = -pi*(u+1) .* sin (pi.*v/2)/2;
  jac(2,2,:) =  pi*(u+1) .* cos (pi.*v/2)/2;
  
  if (nargout == 1)
    varargout{1} = jac;
  elseif (nargout == 2)
    F = ring_polar_map (pts);
    varargout{1} = F;
    varargout{2} = jac;
  end
  
end

