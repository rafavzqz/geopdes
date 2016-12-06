
nel=10;
coeff  = @(x, y,z) ones(size(x));
pmax = 1;

uex     = @(x, y) exp (x) .* sin (y);
h = @(x, y, ind) exp (x) .* sin(y);

for p = pmax:pmax

degree     = [p p];       % Degree of the splines
regularity = [p-1 p-1];       % Regularity of the splines
nsub       = [nel nel];       % Number of subdivisions
nquad      = [p+1 p+1];       % Points for the Gaussian quadrature rule

geometry   = geo_load ('geo_square.txt');
geometry   = geo_load ('geo_ring.txt');
geometry   = geo_load ('geo_roof.txt');

% Questa non va bene per p=1
nrb = nrb4surf([0 0], [1 0], [0 1], [1 2]);

% Questa non ha soluzione esatta
% nrb = nrb4surf([0 0 0], [1 1 0.5], [0 2 0.3], [2 3 3]);
% nrb = nrbdegelev(nrb, [2 2]);
% nrb = nrbmodp (nrb, [0.2 -0.1 -0.6], [6 11]);
% nrb = nrbmodp (nrb, [0.1 -0.3 1.6], [7 10]);
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

rhs = zeros (space.ndof, 1);
rhs_new = rhs;

norm(Stiff_new - Stiff_geopdes,'fro')

[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, [1 2 3 4]);


u_new = zeros (space.ndof, 1); u_geo = zeros (space.ndof, 1);
u_new(drchlt_dofs) = u_drchlt;
u_geo(drchlt_dofs) = u_drchlt;
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs_new(int_dofs) = rhs_new(int_dofs) - Stiff_new(int_dofs, drchlt_dofs)*u_drchlt;
rhs(int_dofs) = rhs(int_dofs) - Stiff_geopdes(int_dofs, drchlt_dofs)*u_drchlt;

u_new(int_dofs) = Stiff_new(int_dofs, int_dofs) \ rhs_new(int_dofs);
u_geo(int_dofs) = Stiff_geopdes(int_dofs, int_dofs) \ rhs(int_dofs);

[eu_new,F] = sp_eval(u_new, space, geometry, [51 51]);
sp_plot_solution(u_new,space,geometry, [51 51]);
figure(2)
[eu_geo,F] = sp_eval(u_geo, space, geometry, [51 51]);
sp_plot_solution(u_geo,space,geometry, [51 51]);
% shading interp

max(abs(u_new(:) - u_geo(:)))
max (abs (eu_new(:) - eu_geo(:)))

figure(3)
sp_plot_solution (u_new - u_geo, space, geometry, [51 51]);
% shading interp

sp_l2_error (space, msh, u_geo-u_new, @(x,y,z) zeros (size(x)))

disp('Errore in norma L2')
sp_l2_error (space, msh, u_new, uex)
sp_l2_error (space, msh, u_geo, uex)

end



