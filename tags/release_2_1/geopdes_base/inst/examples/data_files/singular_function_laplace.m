% SINGULAR_FUNCTION_LAPLACE: auxiliary function to compute the singular solution of the L-shaped domain.

function f = singular_function_laplace (x, y, k)
  [theta, rad] = cart2pol (x, y);
  theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
  f = rad.^(k*2/3) .* sin((k*2/3)*theta);
end

