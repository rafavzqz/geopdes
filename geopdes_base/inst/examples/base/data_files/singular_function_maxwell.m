% SINGULAR_FUNCTION_MAXWELL: auxiliary function to compute the singular solution of the L-shaped domain.

function f = singular_function_maxwell (x, y, k)
  f = zeros(2, size(x,1), size(x,2));
  [theta, rad] = cart2pol (x, y);
  theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
  f(1,:,:) = k*(2/3)*rad.^(k*2/3-1) .* sin((k*2/3-1)*theta);
  f(2,:,:) = k*(2/3)*rad.^(k*2/3-1) .* cos((k*2/3-1)*theta);
end

