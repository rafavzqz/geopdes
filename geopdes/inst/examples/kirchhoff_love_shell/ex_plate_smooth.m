clear all

% Construct geometry structure
base = 10;
p11 = [0 0 1]; p12 = [base 0 1]; p21 = [0 base 1]; p22 = [base base 1];
srf = nrb4surf (p11,p12,p21,p22);

% Physical parameters
E = 1e6;
nu = 0.;
t = (12 * (1 - nu*nu) / E)^(1/3);

% Source and boundary terms
% Distributed load.
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) 4*16*pi^4/(base^4)*sin(2*pi*x/base).*sin(2*pi*y/base);

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));


ux = @(x, y, z) zeros(size(x));
uy = @(x, y, z) zeros(size(x));
uz = @(x, y, z) sin(2*pi*x/base).*sin(2*pi*y/base);

problem_data.uex     = @(x, y, z) cat(1, ...
    reshape (ux (x,y,z), [1, size(x)]), ...
    reshape (uy (x,y,z), [1, size(x)]), ...
    reshape (uz (x,y,z), [1, size(x)]));

problem_data.h     = @(x, y, z, iside) cat(1, ...
    reshape (zeros(size(x)), [1, size(x)]), ...
    reshape (zeros(size(x)), [1, size(x)]), ...
    reshape (zeros(size(x)), [1, size(x)]));

% components of the gradient
ux_dx = @(x, y, z) zeros(size(x));
ux_dy = @(x, y, z) zeros(size(x));
ux_dz = @(x, y, z) zeros(size(x));
uy_dx = @(x, y, z) zeros(size(x));
uy_dy = @(x, y, z) zeros(size(x));
uy_dz = @(x, y, z) zeros(size(x));
uz_dx = @(x, y, z) 2*pi/(base^1)*cos(2*pi*x/base).*sin(2*pi*y/base);
uz_dy = @(x, y, z) 2*pi/(base^1)*sin(2*pi*x/base).*cos(2*pi*y/base);
uz_dz = @(x, y, z) zeros(size(x));


problem_data.graduex = @(x, y, z) cat (1, ...
    reshape (ux_dx(x,y,z), [1, size(x)]), ...
    reshape (uy_dx(x,y,z), [1, size(x)]), ...
    reshape (uz_dx(x,y,z), [1, size(x)]), ...
    reshape (ux_dy(x,y,z), [1, size(x)]), ...
    reshape (uy_dy(x,y,z), [1, size(x)]), ...
    reshape (uz_dy(x,y,z), [1, size(x)]), ...
    reshape (ux_dz(x,y,z), [1, size(x)]), ...
    reshape (uy_dz(x,y,z), [1, size(x)]), ...
    reshape (uz_dz(x,y,z), [1, size(x)]));

% components of the hessian - u_x
ux_dxx = @(x, y, z) zeros(size(x));
ux_dxy = @(x, y, z) zeros(size(x));
ux_dxz = @(x, y, z) zeros(size(x));
ux_dyx = @(x, y, z) zeros(size(x));
ux_dyy = @(x, y, z) zeros(size(x));
ux_dyz = @(x, y, z) zeros(size(x));
ux_dzx = @(x, y, z) zeros(size(x));
ux_dzy = @(x, y, z) zeros(size(x));
ux_dzz = @(x, y, z) zeros(size(x));

% components of the hessian - u_y
uy_dxx = @(x, y, z) zeros(size(x));
uy_dxy = @(x, y, z) zeros(size(x));
uy_dxz = @(x, y, z) zeros(size(x));
uy_dyx = @(x, y, z) zeros(size(x));
uy_dyy = @(x, y, z) zeros(size(x));
uy_dyz = @(x, y, z) zeros(size(x));
uy_dzx = @(x, y, z) zeros(size(x));
uy_dzy = @(x, y, z) zeros(size(x));
uy_dzz = @(x, y, z) zeros(size(x));

% components of the hessian - u_z
uz_dxx = @(x, y, z) -4*pi^2/(base^2)*sin(2*pi*x/base).*sin(2*pi*y/base);
uz_dxy = @(x, y, z) 4*pi^2/(base^2)*cos(2*pi*x/base).*cos(2*pi*y/base);
uz_dxz = @(x, y, z) zeros(size(x));
uz_dyx = @(x, y, z) 4*pi^2/(base^2)*cos(2*pi*x/base).*cos(2*pi*y/base);
uz_dyy = @(x, y, z) -4*pi^2/(base^2)*sin(2*pi*x/base).*sin(2*pi*y/base);
uz_dyz = @(x, y, z) zeros(size(x));
uz_dzx = @(x, y, z) zeros(size(x));
uz_dzy = @(x, y, z) zeros(size(x));
uz_dzz = @(x, y, z) zeros(size(x));

problem_data.hessuex = @(x, y, z) cat (1, ...
    reshape (ux_dxx(x,y,z), [1, size(x)]), ...
    reshape (uy_dxx(x,y,z), [1, size(x)]), ...
    reshape (uz_dxx(x,y,z), [1, size(x)]), ...
    reshape (ux_dyx(x,y,z), [1, size(x)]), ...
    reshape (uy_dyx(x,y,z), [1, size(x)]), ...
    reshape (uz_dyx(x,y,z), [1, size(x)]), ...
    reshape (ux_dzx(x,y,z), [1, size(x)]), ...
    reshape (uy_dzx(x,y,z), [1, size(x)]), ...
    reshape (uz_dzx(x,y,z), [1, size(x)]), ...
    reshape (ux_dxy(x,y,z), [1, size(x)]), ...
    reshape (uy_dxy(x,y,z), [1, size(x)]), ...
    reshape (uz_dxy(x,y,z), [1, size(x)]), ...
    reshape (ux_dyy(x,y,z), [1, size(x)]), ...
    reshape (uy_dyy(x,y,z), [1, size(x)]), ...
    reshape (uz_dyy(x,y,z), [1, size(x)]), ...
    reshape (ux_dzy(x,y,z), [1, size(x)]), ...
    reshape (uy_dzy(x,y,z), [1, size(x)]), ...
    reshape (uz_dzy(x,y,z), [1, size(x)]), ...
    reshape (ux_dxz(x,y,z), [1, size(x)]), ...
    reshape (uy_dxz(x,y,z), [1, size(x)]), ...
    reshape (uz_dxz(x,y,z), [1, size(x)]), ...
    reshape (ux_dyz(x,y,z), [1, size(x)]), ...
    reshape (uy_dyz(x,y,z), [1, size(x)]), ...
    reshape (uz_dyz(x,y,z), [1, size(x)]), ...
    reshape (ux_dzz(x,y,z), [1, size(x)]), ...
    reshape (uy_dzz(x,y,z), [1, size(x)]), ...
    reshape (uz_dzz(x,y,z), [1, size(x)]));


%% Uniform refinement loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deg = 3;
regularity = deg-1;

num_max_iter = 5;

fprintf('Degree: %d \n', deg);

for iter = 1:num_max_iter
    method_data.degree      = [deg deg];        % Degree of the splines
    for icomp=1:3
        method_data.regularity{icomp}  = [regularity regularity];        % Regularity of the splines
    end
    method_data.nsub_coarse = [2^(iter+1) 2^(iter+1)];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
    method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
    method_data.nquad       = [deg+1 deg+1];        % Points for the Gaussian quadrature rule
    method_data.space_type  = 'standard'; % 'simplified' (only children functions) or 'standard' (full basis)
    method_data.truncated   = 0;            % 0: False, 1: True
    
    geometry = geo_load (srf);
    degelev  = max (method_data.degree - (geometry.nurbs.order-1), 0);
    nurbs    = nrbdegelev (geometry.nurbs, degelev);
    [rknots, zeta, nknots] = kntrefine (nurbs.knots, method_data.nsub_coarse-1, nurbs.order-1, [regularity regularity]);
    
    nurbs = nrbkntins (nurbs, nknots);
    geometry = geo_load (nurbs);
    
    problem_data.geo_name = nurbs;
    
    % Construct msh structure
    rule     = msh_gauss_nodes (method_data.nquad);
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh      = msh_cartesian (zeta, qn, qw, geometry,'der2', true);
    
    problem_data.E_coeff = @(x, y, z) E * ones(size(x));
    problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
    problem_data.t_coeff = @(x, y, z) t * ones(size(x));
    
    % Construct the reference space.
    % sp_scalar = sp_nurbs (problem_data.geo_name, msh);
    sp_scalar = sp_bspline (rknots, method_data.degree, msh);
    scalar_spaces = repmat ({sp_scalar}, 1, msh.rdim);
    space = sp_vector (scalar_spaces, msh);
    
    
    ndof = space.ndof;
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',iter);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SOLVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('SOLVE:')
    % Assemblying stiffness matrix ...
    K = op_KL_shells_tp (space, space, msh, problem_data.E_coeff, problem_data.nu_coeff, problem_data.t_coeff);
    % ... and right-hand-side
    rhs = op_f_v_tp (space, msh, problem_data.f);
    
    u = zeros (ndof, 1);
    
    % Apply pinned boundary conditions
    drchlt_dofs = space.boundary(1).dofs(space.boundary(1).comp_dofs{3});
    drchlt_dofs = union(drchlt_dofs,space.boundary(2).dofs(space.boundary(2).comp_dofs{3}));
    drchlt_dofs = union(drchlt_dofs,space.boundary(3).dofs(space.boundary(3).comp_dofs{3}));
    drchlt_dofs = union(drchlt_dofs,space.boundary(4).dofs(space.boundary(4).comp_dofs{3}));
    
    drchlt_dofs = union(drchlt_dofs,space.boundary(1).dofs(space.boundary(1).comp_dofs{2}));
    drchlt_dofs = union(drchlt_dofs,space.boundary(2).dofs(space.boundary(2).comp_dofs{2}));
    drchlt_dofs = union(drchlt_dofs,space.boundary(3).dofs(space.boundary(3).comp_dofs{1}));
    drchlt_dofs = union(drchlt_dofs,space.boundary(4).dofs(space.boundary(4).comp_dofs{1}));
    
    int_dofs = setdiff (1:ndof, drchlt_dofs);
    rhs(int_dofs) = rhs(int_dofs) - K(int_dofs, drchlt_dofs)*u(drchlt_dofs);
    
    % Solve the linear system
    u(int_dofs) = K(int_dofs, int_dofs) \ rhs(int_dofs);
    fprintf('Number of elements: %d. Total DOFs: %d \n', msh.nel, space.ndof);
    
    energy_norm(iter) = u'*K*u;
    fprintf('Energy norm: %d. \n', energy_norm(iter));
    
    
    if (isfield (problem_data, 'hessuex'))
        msh_eval = msh_precompute(msh);
        space_eval = sp_precompute(space, msh_eval, 'gradient', true, 'hessian', true);
        h(iter) = min(msh_eval.element_size);
        [~, ~, err_l2(iter), err_h2s(iter), err_h1s(iter), ~, ~, ~, ~, ~] = sp_h2_error (space_eval, msh_eval, u, ...
            problem_data.uex, problem_data.graduex, problem_data.hessuex);
    end
    
    fprintf('Error in L2 norm = %g\n', err_l2(iter));
    fprintf('Error in H1 semi norm = %g\n', err_h1s(iter));
    fprintf('Error in H2 semi norm = %g\n', err_h2s(iter));
    
    
end

figure()
factor = 5e-3;
loglog(h, err_h2s, '-x', h, factor * h.^(deg-1)); legend('error H^2', strcat('slope h^',num2str(deg-1)), 'location', 'best')
figure()
loglog(h, err_h1s, '-x', h, factor * h.^(deg)); legend('error H^1', strcat('slope h^',num2str(deg)), 'location', 'best')
figure()
loglog(h, err_l2, '-x', h, factor * h.^(deg+1)); legend('error L^2', strcat('slope h^',num2str(deg+1)), 'location', 'best')

disp('END:')






