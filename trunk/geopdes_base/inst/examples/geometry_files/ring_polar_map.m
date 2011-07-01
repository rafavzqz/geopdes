% RING_POLAR_MAP: parameterization of a ring using polar coordinates.

function F = ring_polar_map (pts)
  
  if (iscell (pts))
    u = reshape (repmat (pts{1}(:), 1, numel (pts{2})), 1, []);
    v = reshape (repmat (pts{2}(:)', numel (pts{1}), 1), 1, []);
  else
    u = pts(1,:); v = pts(2,:);
  end
  
  F(1,:) = (u+1) .* cos (pi.*v/2);
  F(2,:) = (u+1) .* sin (pi.*v/2);
end

