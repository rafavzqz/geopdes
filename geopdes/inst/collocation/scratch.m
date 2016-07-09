clear 

%% ================================ problem and choice of discretization  ================================

% problem_case = 1; % square 2D
problem_case = 2; % ring 2D

switch problem_case
    
    case 1    
        % Physical domain, defined as NURBS map given in a text file
        problem_data.geo_name = 'geo_square.txt';        
        % initial import of geometry
        geom_no_ref  = geo_load (problem_data.geo_name);
        problem_data.D = geom_no_ref.rdim;
        
        % Type of boundary conditions for each side of the domain
        problem_data.nmnn_sides   = [];
        problem_data.drchlt_sides = [1 2 3 4];        
        
        % Physical parameters
        problem_data.c_diff  = @(x, y) ones(size(x));
        
        % Source and boundary terms
        problem_data.f = @(x, y) 2*sin(pi*y) + pi^2*x.*(1-x).*sin(pi*y);
        problem_data.g = @(x, y, ind) zeros(size(x));
        problem_data.h = @(x, y, ind) zeros(size(x));
        
        % Exact solution (optional)
        problem_data.uex     = @(x, y) x.*(1-x) .* sin(pi*y);
        problem_data.graduex = @(x, y) cat (1, ...
            reshape ((1-2*x).*sin(pi*y), [1, size(x)]), ...
            reshape (pi*x.*(1-x).*cos(pi*y), [1, size(x)]));
        
    case 2
        % Physical domain, defined as NURBS map given in a text file
        problem_data.geo_name = 'geo_ring.txt';        
        % initial import of geometry
        geom_no_ref  = geo_load (problem_data.geo_name);
        geom_no_ref.nurbs = nrbdegelev(geom_no_ref.nurbs,[1 0]); % this way I get degree two in each direction        
        problem_data.D = geom_no_ref.rdim;

        % Type of boundary conditions for each side of the domain
        problem_data.nmnn_sides   = [];
        problem_data.drchlt_sides = [1 2 3 4];        

        % Physical parameters
        problem_data.c_diff  = @(x, y) ones(size(x));
        
        % Source and boundary terms
        problem_data.f = @(x,y) 2*x.*(22.*x.^2.*y.^2+21.*y.^4-45.*y.^2+x.^4-5.*x.^2+4);
        problem_data.g = @(x, y, ind) zeros(size(x));
        problem_data.h = @(x, y, ind) zeros(size(x));
        
        % Exact solution (optional)
        problem_data.uex     = @(x, y) -(x.^2+y.^2-1).*(x.^2+y.^2.-4).*x.*y.^2;
        problem_data.graduex = @(x, y) cat (1,  reshape (-2*(x.*y).^2.*((x.^2+y.^2-1)+(x.^2+y.^2-4)) - (x.^2+y.^2-1).*(x.^2+y.^2-4).*y.^2, [1, size(x)]), ...
                                        reshape ( -2*x.*y.^3.*((x.^2+y.^2-1)+(x.^2+y.^2-4)) -   2*x.*y.*(x.^2+y.^2-1).*(x.^2+y.^2-4), [1, size(x)]));        
                
end
                   
                   
                   
% discretization parameters (p and h)

method_data.degree     = [3 3];       % Degree of the splines, obtained by k-refinement of geometry
method_data.n_sub      = [16 16];     % will divide each subinterval of the original knot span in n_sub many subinterval
                                      % i.e., we add nsub-1 knots in each interval of the original knotline.
                                      % note that if the original geometry has even number of subintervals, all possible
                                      % refinements with this strategy will have even subintervals

%% ================================  geometry ================================


% make sure degree is identical in each direction
if length(unique(geom_no_ref.nurbs.order))>1
    error('different degrees in different dim')
end

% compute how many degrees to elev. If 0 or less, something is wrong
degelev  = max (method_data.degree - (geom_no_ref.nurbs.order-1), 0);
if max(degelev)<=0,
    error(strcat('max(degelev)=',num2str(max(degelev)),', stop! You need to fix the user set simulation degree'))
end
% otherwise, do degree elev
nurbs_ptemp    = nrbdegelev (geom_no_ref.nurbs, degelev);

% next h-ref (so the overall procedure is k-ref). 
new_knots = cell(1,problem_data.D);

for d=1:problem_data.D
        
    knotline = unique(nurbs_ptemp.knots{d});
    nb_sub0 = length(knotline)-1;
    NN = method_data.n_sub(d);
    new_knots{d} = [];

    for i = 1:nb_sub0
        tmp_new = linspace(knotline(i),knotline(i+1),NN+1);
        tmp_new(1)=[]; tmp_new(end)=[];
        new_knots{d}(end+1:end+method_data.n_sub(d)-1) = tmp_new;
    end
    
end

% this is the final domain to use for computations
nurbs = nrbkntins( nurbs_ptemp, new_knots);
geo_refined = geo_load (nurbs);

% these are the knot lines (with repetitions)
knots=geo_refined.nurbs.knots;


%% ================================ Collocation ================================


% now we need to generate the collocation points. We define a sequence of coll_pts and create
% a fictitious knot line so that we have 1 point per element. Then we build a mesh object
% over this knot line and the collocation points. Then Geopdes treats them as quadrature points, 
% but provides also a mean of evaluating the basis functions in those points, we go for it


%pts_case = 1; % equispaced
pts_case = 2; % greville -----------------------------------------------------------------> with greville I have coll-pts = DoFs so in 
                                                                                        % principle I do not need to use least squares. 
                                                                                        % However, I keep using it for the moment because 
                                                                                        % I want to be general in the choice of points
                                                                                        % and also the Dir BC are treated eliminating 
                                                                                        % the DoFs (matrix columns) but not rows
switch pts_case
    case 1
        cpt_1=linspace(0,1,40); cpt_1(1)=[]; cpt_1(end)=[];
        cpt_2=linspace(0,1,40); cpt_2(1)=[]; cpt_2(end)=[];
    case 2
        %  aveknt(knot_line,k)  returns the averages of successive  k-1  knots. To comply
        % with the definition of greville as g_j = ( zeta_{j+1} + ... + zeta_{j+p})/p
        % for an open knot line with n+p+1 and degree p, we then pass as second argument p+1
        cpt_1 = aveknt(knots{1},method_data.degree(1)+1); 
        cpt_2 = aveknt(knots{2},method_data.degree(2)+1); 
end

coll_pts={cpt_1,cpt_2};


% for each dim, we generate a fictitious knot line, so that there is one coll point per element.
% we define this by taking as knot line [0 midpoint(coll-pt1, coll-pt2) midpoint(coll-pt2, coll-pt3) ... 1]

for d = 1:problem_data.D
    coll_pts{d} = coll_pts{d}(:)';
    if numel(coll_pts{d}) > 1
        brk{d} = [knots{d}(1), coll_pts{d}(1:end-1) + diff(coll_pts{d})/2, knots{d}(end)];
    else
        brk{d} = [knots{d}(1) knots{d}(end)];
    end
end

% here's the mesh with one coll pt per element
coll_msh = msh_cartesian (brk, coll_pts, [], geo_refined, 'boundary', true,'der2',true);

% now we generate the spline basis function set
space = sp_nurbs( geo_refined.nurbs, coll_msh );


% now we evaluate each spline in the collocation points and collect everything into the design matrix, i.e. a
% matrix whose n-th row is the equation collocated in the n-th point. Remember that we have built coll_msh
% to have 1 pt per element, so the number of coll pts is exactly the number of elements of coll_msh

tot_nb_coll_pts = coll_msh.nel;
nb_dofs = space.ndof;
A = zeros(tot_nb_coll_pts,nb_dofs);

% evaluate parameterization 
coll_mesh_eval = msh_precompute(coll_msh);
% evaluate basis fun 
sp_evals  = sp_precompute(space, coll_mesh_eval, 'gradient', false, 'laplacian', true);

% sp_evals.shape_function_laplacian contains evaluations stored in a matrix of size
% [msh.nqn x nsh_max x msh.nel double]      
% i.e. for each node of an element (msh.nqn) only nonzero functions (nsh_max) 
% the nonzero functions are included in the matrix connectivity, whose size is (nsh_max x msh.nel vector)         
% so we need to unpack this information. To do this, we loop over the number of elements (i.e., over coll pts for this case),
% detect the list of nonzero functions on that element and store the evaluations in A

for iel = 1:coll_msh.nel %remember that although this formally spans elements, in practice it spans coll pts too
    
    %the nonzero funs on this element
    list_fun = sp_evals.connectivity(:,iel);
    
    % since we have only 1 pt per element, shape_function_laplacian has in practice size
    % [1 x nonzerofun x nb_mesh_el] 
    A(iel,list_fun) = -sp_evals.shape_function_laplacians(1,:,iel); % sidenote ------------------------------> not this easy if instead we are solving -div{a grad u}=f
    
end


% generate rhs. Recover coordinates of collocation points and evaluate forcing 
x={};
for idir = 1:coll_mesh_eval.rdim
    x{idir} = reshape (coll_mesh_eval.geo_map(idir,:,:), coll_mesh_eval.nqn*coll_msh.nel, 1);
end
rhs  = problem_data.f(x{1},x{2});


% remove DoF of boundary condition. Sidenote ----------------------------------------------------------> homogeneous Dir hardcoded here

nb_boundaries = length(sp_evals.boundary);
boundary_dofs=[];
for bb=1:nb_boundaries
    boundary_dofs = union(boundary_dofs,sp_evals.boundary(bb).dofs);
end
A(:,boundary_dofs)=[];
internal_dofs = setdiff(1:sp_evals.ndof,boundary_dofs);

% solve Ax=rhs sidenote -------------------------------------------------------------------------------> with least squares for now A'*A 
AA = A'*A;
rr = A'*rhs;

sol = AA\rr;

% finally, put together boundary and internal_dofs
u_coll = zeros(sp_evals.ndof,1);

u_coll(internal_dofs)=sol;



%% ================================ Galerkin for comparison ================================

% fix quad points, for good this time, and build mesh
nquad      = [5 5];     % Points for the Gaussian quadrature rule
[gal_qn, gal_qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (nquad));
gal_msh      = msh_cartesian (knots, gal_qn, gal_qw, geo_refined);
  
% Even if we want to use the same space as before, we need to rebuild the space object because now the mesh is
% different
gal_space = sp_nurbs(geo_refined.nurbs, gal_msh);
% gal_space    = sp_bspline (knots, method_data.degree, gal_msh);


% Assemble the matrices
gal_stiff_mat = op_gradu_gradv_tp (gal_space, gal_space, gal_msh, problem_data.c_diff);
gal_rhs       = op_f_v_tp (gal_space, gal_msh, problem_data.f);

% Apply Dirichlet boundary conditions --------------------------------------------------------------------> homogenous dir hardcoded
u_gal = zeros (gal_space.ndof, 1);

gal_nb_boundaries = length(gal_space.boundary);
gal_boundary_dofs=[];
for bb=1:gal_nb_boundaries
    gal_boundary_dofs = union(gal_boundary_dofs,gal_space.boundary(bb).dofs);
end
gal_int_dofs = setdiff (1:gal_space.ndof, gal_boundary_dofs);


% Solve the linear system
u_gal(gal_int_dofs) = gal_stiff_mat(gal_int_dofs, gal_int_dofs) \ gal_rhs(gal_int_dofs);




%% ============================== the comparison =======================================


% plot of solution
plot_pts = {linspace(0, 1, 40), linspace(0, 1, 40)};
figure
[eu, F] = sp_eval (u_coll, space, geo_refined, plot_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,3,1)
surf (X, Y, eu)
title ('collocation solution'), axis tight
subplot (1,3,2)
[eu, F] = sp_eval (u_gal, gal_space, geo_refined, plot_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('galerkin solution'), axis tight
subplot (1,3,3)
surf (X, Y, problem_data.uex (X,Y))
title ('Exact solution'), axis tight



% compute errors of coll and gal. Create yet another mesh, which provides quad points for the error
quad_err_nquad      = [9 9];     % Points for the Gaussian quadrature rule
[quad_err_qn, quad_err_qw] = msh_set_quad_nodes (knots, msh_gauss_nodes(quad_err_nquad));
quad_err_msh      = msh_cartesian (knots, quad_err_qn, quad_err_qw, geo_refined);

% I also build another space, to evaluate the basis functions in the right quadrature points for error
quad_err_space = sp_nurbs( geo_refined.nurbs, quad_err_msh ); 


% Display errors of the computed solution in the L2 and H1 norm
format long
[error_h1_coll, error_l2_coll] = sp_h1_error (quad_err_space, quad_err_msh, u_coll, problem_data.uex, problem_data.graduex);
[error_h1_gal, error_l2_gal] = sp_h1_error (quad_err_space, quad_err_msh, u_gal, problem_data.uex, problem_data.graduex);

[error_l2_gal error_h1_gal error_l2_coll error_h1_coll]



%% ======= a small table of numerical errors obtained by repeatedly running the code with different h, to check convergence order ======

% settings: ring problem 2D, greville abscissae, spline degree 3, n_sub = [4 8 16 32 64]
% [error_l2_gal error_h1_gal error_l2_coll error_h1_coll]
err= [ 0.006607242321582   0.098667163584905   0.221683362376558   0.826978587598251;
       0.000321871231424   0.012273598137088   0.070945374842625   0.263695190251720;
       0.000019497341911   0.001588947267174   0.019005325968532   0.071420226306507;
       0.000001224267705   0.000203182231291   0.004839925626712   0.018288070038202;
       0.000000077011984   0.000025712824749   0.001215783435746   0.004602652148439;
];

h = 1./[4 8 16 32 64];    

figure
loglog(h,err(:,1),'-ob','DisplayName','L^2 Galerkin')
hold on
loglog(h,h.^4,'--b','DisplayName','h^4')
loglog(h,err(:,2),'-or','DisplayName','H^1 Galerkin')
loglog(h,h.^3,'--r','DisplayName','h^3')
loglog(h,err(:,3),'-ok','DisplayName','L^2 Coll Greville')
loglog(h,h.^2,'--k','DisplayName','h^2')
loglog(h,err(:,4),'-x','Color',[0 0.8 0],'DisplayName','H^1 Coll Greville')
grid on
legend show
set(legend,'Location','SouthEast')
