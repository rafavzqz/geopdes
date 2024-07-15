% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';
    
% Physical parameters
lambda = (1/(4*sqrt(2)*pi))^2;
problem_data.lambda  = @(x, y) lambda* ones(size(x));

% Periodic directions and boundary conditions
problem_data.periodic_directions = [];

% Time and time step size
Time_max = .08;
dt = 1e-3;
time_step_save = linspace(dt,Time_max,20);
problem_data.time = 0;
problem_data.Time_max = Time_max;


% Penalty parameters
problem_data.Cpen_nitsche = 1e4 * lambda; % Nitsche's method parameter
problem_data.Cpen_projection = 1000;      % parameter of the penalized L2 projection (see initial conditions)

% 2) INITIAL CONDITIONS
% mean = 0.0;
% var = 0.2;
% ic_fun = @(x, y) mean + (rand(size(x))*2-1)*var; % Random initial condition
ic_fun = @(x, y) 0.1 * cos(2*pi*x) .* cos(2*pi*y); % Condition as in Gomes, Reali, Sangalli, JCP (2014).
problem_data.fun_u = ic_fun;
% problem_data.fun_udot = [];

% 3) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
deg = 2; nel = 10;
method_data.degree     = [deg deg];    % Degree of the splines
method_data.regularity = [deg-1 deg-1];    % Regularity of the splines
method_data.nsub       = [nel nel];  % Number of subdivisions
method_data.nquad      = [deg+1 deg+1];    % Points for the Gaussian quadrature rule

% time integration parameters
method_data.rho_inf_gen_alpha = 0.5; % Parameter for generalized-alpha method
method_data.dt = dt;                 % Time step size



% 4) CALL TO THE SOLVER
[geometry, msh, space, results] = solve_cahn_hilliard (problem_data, method_data, time_step_save);

% 5) POST-PROCESSING        
vtk_pts = {linspace(0, 1, nel*4), linspace(0, 1, nel*4)};
folder_name = strcat('results_CH_p',num2str(deg),'_nel',num2str(nel),'_lambda',num2str(lambda));
status = mkdir(folder_name);

% 5.1) EXPORT TIME  
filename = strcat( folder_name,'/filenum_to_time.mat');
time_steps = results.time;
save(filename, 'time_steps');


% 5.2) EXPORT TO PARAVIEW    
for step = 1:length(results.time)
  output_file = strcat( folder_name,'/Square_cahn_hilliard_', num2str(step) );
  fprintf ('The result is saved in the file %s \n \n', output_file);
  sp_to_vtk (results.u(:,step), space, geometry, vtk_pts, output_file, {'u','grad_u'}, {'value','gradient'})
end
    
% 5.3) PLOT LAST RESULT
[eu, F] = sp_eval (results.u(:,end), space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
colorbar
view(0,90)
shading interp
axis equal tight