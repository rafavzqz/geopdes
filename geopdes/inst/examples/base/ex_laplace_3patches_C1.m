% PHYSICAL DATA OF THE PROBLEM
clear problem_data
problem_data.geo_name = 'geo_3patch_ASG1.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.weak_drchlt_sides = [1 2 3 4 5 6];

problem_data.c_diff  = @(x, y) ones(size(x));
        
% Source and boundary terms
C = 20; 
P = [17/3+0.15, 2+0.075];
alpha = 4/6;
problem_data.f = @(x,y) -2*alpha*((x-P(1)).^2+(y-P(2)).^2).^(alpha-2).*...
                       (((x-P(1)).^2+(y-P(2)).^2)+2*(alpha-1)*(x-P(1)).^2)...
                        -2*alpha*((x-P(1)).^2+(y-P(2)).^2).^(alpha-2).*...
                       (((x-P(1)).^2+(y-P(2)).^2)+2*(alpha-1)*(y-P(2)).^2);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) ((x-P(1)).^2+(y-P(2)).^2).^alpha;

% Exact solution (optional)
problem_data.uex =@(x,y) ((x-P(1)).^2+(y-P(2)).^2).^alpha;
problem_data.graduex = @(x,y) 2*alpha*cat (1, ...
            reshape (((x-P(1)).^2+(y-P(2)).^2).^(alpha-1).*(x-P(1)), [1, size(x)]), ...
            reshape (((x-P(1)).^2+(y-P(2)).^2).^(alpha-1).*(y-P(2)), [1, size(x)]));

% DISCRETIZATION PARAMETERS
clear method_data
deg = 4;
method_data.degree      = [deg deg];     % Degree of the splines
method_data.regularity  = [deg-2 deg-2]; % Regularity of the splines
method_data.nsub        = [8 8];         % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad       = [deg+1 deg+1]; % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);

[geometry, msh, space, u] = mp_solve_laplace_C1 (problem_data, method_data);
err_h1 = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)

npts = [81 81];
sp_plot_solution (u, space, geometry, npts); shading interp
