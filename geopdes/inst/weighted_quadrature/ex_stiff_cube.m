
nel=20;
coeff  = @(x, y, z) ones(size(x));
pmax = 2;

for p = pmax:pmax

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

tempo = tic;
Stiff_new = Stiff_WQ(msh, space, geometry, coeff);
tempo_nostro = toc(tempo)

tempo = tic;
Stiff_geopdes = op_gradu_gradv_tp(space, space, msh, coeff);
tempo_geopdes = toc(tempo)

norm(Stiff_new - Stiff_geopdes,'fro')

end



