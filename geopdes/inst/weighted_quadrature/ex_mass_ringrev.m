deg = 3;
nel = 16;
coeff  = @(x, y, z) ones(size(x));

% Function to project
f = @(x,y,z) cos(x) .* sqrt(y) .* log(y.*z+1);

for ii = 1:4
  nel = 2^ii;
  degree     = [deg deg deg];       % Degree of the splines
  regularity = [deg-1 deg-1 deg-1];       % Regularity of the splines
  nsub       = [nel nel nel];       % Number of subdivisions
  nquad      = [deg+1 deg+1 deg+1];       % Points for the Gaussian quadrature rule

%%%%%%%%%%% creo la geometria %%%%%%%%%%%%%%
  geometry = geo_load ('geo_ring.txt');
  vol = nrbrevolve (geometry.nurbs, [-1 -1 -1], [0 -1 0], pi/4);
  geometry = geo_load (vol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
  rule     = msh_gauss_nodes (nquad);
  [qn, qw] = msh_set_quad_nodes (zeta, rule);
  msh      = msh_cartesian (zeta, qn, qw, geometry);
  space    = sp_bspline (knots, degree, msh);

  Mass_new = Mass_WQ (msh,space,geometry, coeff);
  Mass_geopdes = op_u_v_tp (space, space, msh, coeff);
  
  rhs = op_f_v_tp (space, msh, f);
  err_mat(ii) = norm(Mass_geopdes - Mass_new,'fro')/norm(Mass_geopdes,'fro');
  
  u_new = Mass_new \ rhs;
  u_geopdes = Mass_geopdes \ rhs;
  
  err_new(ii) = sp_l2_error (space, msh, u_new, f);
  err_geopdes(ii) = sp_l2_error (space, msh, u_geopdes, f);
  
  ndof(ii) = space.ndof;
end

