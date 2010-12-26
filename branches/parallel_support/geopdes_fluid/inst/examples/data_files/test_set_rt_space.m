% TEST_SET_RT_SPACE: auxiliary script to force the use of RT discrete space.

fun_space    = 'sp_bspline_rt_2d_phys';
der2         = true;
if (~isempty (drchlt_sides) && any (degree ~= 3))
  degree     = [3 3];
  regularity = min (regularity, 2); 
  fprintf ('for RT elements with Dirichlet conditions only degree 3 may be used\n')
end
