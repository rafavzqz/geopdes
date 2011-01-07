geometry = geo_load (eye (4));
breaks = {(linspace (0, 1, 11)), (linspace (0, 1, 11))};
knots = kntbrkdegreg (breaks, [4 4], [3 3]);
[qn, qw] = msh_set_quad_nodes (breaks, msh_gauss_nodes ([9 9]));
msh = msh_2d_tensor_product (breaks, qn, qw); 
msh = msh_push_forward_2d (msh, geometry, 'der2', true);
sp = sp_bspline_2d_phys (knots, [4 4], msh, 'hessian', true);

# f = @(x, y) (x.^3 + 3*y.^2 +1);
# dfdx = @(x, y) (3 * x.^2);
# dfdy = @(x, y) (6*y);

f = @(x, y) (x.^2.*(x-1).^2.*y.^2.*(y-1).^2);
dfdx = @(x, y) (2*(x.*(x-1).^2+x.^2.*(x-1)).*y.^2.*(y-1).^2);
dfdy = @(x, y) (2*(y.*(y-1).^2+y.^2.*(y-1)).*x.^2.*(x-1).^2);

df = @(x, y) cat(1, ...
                reshape (dfdx (x,y), [1, size(x)]), ...
                reshape (dfdy (x,y), [1, size(x)]));

# d2fdx2 = @(x, y) (6*x);
# d2fdxy = @(x, y) (zeros (size (x)));
# d2fdy2 = @(x, y) (6 * ones (size (y)));

d2fdx2 = @(x, y) (2*((x-1).^2+4*x.*(x-1)+x.^2).*(y.^2).*(y-1).^2);
d2fdxy = @(x, y) (4*(x.*(x-1).^2+x.^2.*(x-1)).*(y.*(y-1).^2+y.^2.*(y-1)));
d2fdy2 = @(x, y) (2*((y-1).^2+4*y.*(y-1)+y.^2).*(x.^2).*(x-1).^2);

[x, y] = deal (squeeze (msh.geo_map(1,:,:)), squeeze (msh.geo_map(2,:,:)));
M = op_u_v (sp, sp, msh, ones (size (x)));
b = op_f_v (sp, msh, f (x, y));

u = M \ b;

% [eu, F] = sp_eval_2d (u, sp, geometry, {linspace(0,1,11)', linspace(0,1,11)'});
% msh_to_vtk_2d (F, eu, 'pippo.vts', 'pippo');

fprintf ('errl2 = %g\n', sp_l2_error (sp, msh, u, f))
fprintf ('errh1 = %g\n',sp_h1_error (sp, msh, u, f, df))

weights = reshape (u(sp.connectivity), 1, sp.nsh_max, msh.nel);
weights = repmat (weights, [msh.nqn, 1, 1]);

dudx = squeeze (sum (weights .* squeeze (sp.shape_function_gradients(1,:,:,:)), 2)); 

d2udx2 = squeeze (sum (weights .* squeeze (sp.shape_function_hessians(1,1,:,:,:)), 2)); 
fprintf ('errl2_d2udx2 = %g\n',norm (d2udx2 - squeeze (d2fdx2 (x, y))))

d2udxy = squeeze (sum (weights .* squeeze (sp.shape_function_hessians(1,2,:,:,:)), 2)); 
fprintf ('errl2_d2udxy = %g\n',norm (d2udxy - squeeze (d2fdxy (x, y))))

d2udxy = squeeze (sum (weights .* squeeze (sp.shape_function_hessians(2,1,:,:,:)), 2)); 
fprintf ('errl2_d2udxy = %g\n',norm (d2udxy - squeeze (d2fdxy (x, y))))

d2udy2 = squeeze (sum (weights .* squeeze (sp.shape_function_hessians(2,2,:,:,:)), 2)); 
fprintf ('errl2_d2udy2 = %g\n',norm (d2udy2 - squeeze (d2fdy2 (x, y))))
