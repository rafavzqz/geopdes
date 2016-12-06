
nel=160;
coeff  = @(x, y,z) ones(size(x));
pmax = 1;

for p = pmax:pmax

degree     = [p p];       % Degree of the splines
regularity = [p-1 p-1];       % Regularity of the splines
nsub       = [nel nel];       % Number of subdivisions
nquad      = [p+1 p+1];       % Points for the Gaussian quadrature rule

geometry   = geo_load ('geo_square.txt');
geometry   = geo_load ('geo_ring.txt');
geometry   = geo_load ('geo_roof.txt');

% nrb = geometry.nurbs;
% nrb = nrbtform (nrb, vecrot (pi/5, [-1 2 -1]));
nrb = nrb4surf([0 0 0], [1 1 0.5], [0 2 0.3], [2 3 3]);
nrb = nrbdegelev(nrb, [2 2]);
nrb = nrbmodp (nrb, [0.2 -0.1 -0.6], [6 11]);
nrb = nrbmodp (nrb, [0.1 -0.3 1.6], [7 10]);
geometry = geo_load (nrb);

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
Mass = op_u_v_tp(space, space, msh, coeff);
tempo_geopdes = toc(tempo)

rhs = op_f_v_tp (space, msh, coeff);

norm(Stiff_new - Stiff_geopdes,'fro')

u_new = zeros (space.ndof, 1); u_geo = zeros (space.ndof, 1);
u_new(1:space.ndof-1) = Stiff_new(1:end-1, 1:end-1) \ rhs(1:end-1);
u_geo(1:space.ndof-1) = Stiff_geopdes(1:end-1,1:end-1) \ rhs(1:end-1);

[eu_new,F] = sp_eval(u_new, space, geometry, [51 51]);
% surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)), eu_new); shading interp
% figure(2)
[eu_geo,F] = sp_eval(u_geo, space, geometry, [51 51]); 
% shading interp
% surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)), eu_geo);  shading interp

max(abs(u_new(:) - u_geo(:)))
max (abs (eu_new(:) - eu_geo(:)))

% figure(3)
[eu_diff,F] = sp_eval(u_geo-u_new, space, geometry, [51 51]); 
% shading interp
% surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)), eu_diff);  shading interp

sp_l2_error (space, msh, u_geo-u_new, @(x,y,z) zeros (size(x)))

end



