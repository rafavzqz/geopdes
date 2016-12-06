
p=4;
nel=10;

problem_data.geo_name = 'geo_cube.txt';
% Type of boundary conditions for each side of the domain
coeff  = @(x, y, z) ones(size(x));

degree     = [p p p];       % Degree of the splines
regularity = [p-1 p-1 p-1];       % Regularity of the splines
nsub       = [nel nel nel];       % Number of subdivisions
nquad      = [p+1 p+1 p+1];       % Points for the Gaussian quadrature rule
geometry   = geo_load ('geo_cube.txt');
[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);
space    = sp_bspline (knots, degree, msh);

Mass_new = Mass_WQ (msh,space,geometry, coeff);

%Mass_geopdes = op_u_v_tp_modif(space, space, msh, coeff);
%norm(Mass_geopdes - Mass_new,'fro')/norm(Mass_geopdes,'fro')


