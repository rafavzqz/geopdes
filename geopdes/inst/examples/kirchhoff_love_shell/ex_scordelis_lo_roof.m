clear all

E = 4.32e8;
nu = 0.0;
t = 0.25;
% Construct geometry structure
radius = 25.0;
theta = 2.0 / 9.0 * pi;
phi = pi/2 - theta;

pnts = zeros(4,3,2);
w2 = cos(theta);
pnts(:,:,1) = [ -radius*cos(phi)  0.0  radius*cos(phi);
    0.0  0.0  0.0;
    radius*sin(phi)  radius/(cos(theta))*w2  radius*sin(phi);
    1.0  w2  1.0];
pnts(:,:,2) = [ -radius*cos(phi)  0.0 radius*cos(phi);
    2*radius  2*radius*w2  2*radius;
    radius*sin(phi)  radius/(cos(theta))*w2  radius*sin(phi);
    1.0  w2  1.0];
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];

srf = nrbmak(pnts,knots);

problem_data.thickness = t;

% Source and boundary terms
% Distributed load.
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) -90*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

num_max_iter = 5;

plot_flag = false;
%% Uniform refinement loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deg = 3;
regularity = deg-1;

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
    
    % Apply rigid diaphragm boundary conditions    
    drchlt_dofs = space.boundary(3).dofs(space.boundary(3).comp_dofs{1});
    drchlt_dofs = union(drchlt_dofs,space.boundary(3).dofs(space.boundary(3).comp_dofs{3}));
    drchlt_dofs = union(drchlt_dofs,space.boundary(4).dofs(space.boundary(4).comp_dofs{1}));
    drchlt_dofs = union(drchlt_dofs,space.boundary(4).dofs(space.boundary(4).comp_dofs{3}));
    
    int_dofs = setdiff (1:ndof, drchlt_dofs);
    rhs(int_dofs) = rhs(int_dofs) - K(int_dofs, drchlt_dofs)*u(drchlt_dofs);
    
    % Solve the linear system
    u(int_dofs) = K(int_dofs, int_dofs) \ rhs(int_dofs);
    fprintf('Number of elements: %d. Total DOFs: %d \n', msh.nel, space.ndof);
    
    energy_norm(iter) = u'*K*u;
    fprintf('Energy norm: %d. \n', energy_norm(iter));
    
    % To solve to VTK
    if(plot_flag)
        output_file = strcat('scordelisLo_p',num2str(deg),'_GR',num2str(iter));
        vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
        fprintf ('The result is saved in the file %s \n \n', output_file);
        sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
    end
    
    pts{1} = 0.0; pts{2} = 0.5;
    [eu,~] = sp_eval(u,space, geometry,pts);
    reference_displacement = -0.3005925;
    fprintf('Numerical displacement at edge: %d. \n', eu(3));
    fprintf('Reference displacement at edge: %d. \n', reference_displacement);
    
end

disp('END:')


