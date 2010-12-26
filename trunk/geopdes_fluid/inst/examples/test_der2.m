geometry = geo_load (eye(4));
breaks = {(linspace (0, 1, 11)), (linspace (0, 1, 11))};
knots = kntbrkdegreg (breaks, [3 3], [2 2]);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes ([5 5]));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', true);
sp = sp_bspline_2d_phys (knots, [3 3], msh, 'hessian', true);

f = @(x, y) (x.^2 + 2*y);
d2fdx2 = @(x, y) (2 * ones (size (x)));

[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
M = op_u_v (sp, sp, msh, ones (size (x)));
b = op_f_v (sp, msh, f (x, y));

u = M \ b;

[eu, F] = sp_eval_2d (u, sp, geometry, {linspace(0,1,11)', linspace(0,1,11)'});
msh_to_vtk_2d (F, eu, 'pippo.vts', 'pippo');


weights = reshape (u(sp.connectivity), 1, sp.nsh_max, msh.nel);
weights = repmat (weights, [msh.nqn, 1, 1]);

d2udx2 = squeeze (sum (weights .* squeeze (sp.shape_function_hessians(1,1,:,:,:)), 2)); 
norm (d2udx2 - squeeze (d2fdx2 (x, y)))
