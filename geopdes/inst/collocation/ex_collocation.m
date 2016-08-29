% EX_COLLOCATION: 
% This program implements stochastic collocation for the Poisson equation 
% with homogeneous Dirichlet boundary conditions on a 2D/3D domain.
% By Rafael Vazquez and Lorenzo Tamellini
% 
%
% Copyright (C) 2016 Lorenzo Tamellini, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.



%% ================================ problem and choice of discretization  ================================

% problem_case = 1; % square 2D
problem_case = 2; % ring 2D
% problem_case = 3; % ring 3D


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
%         geom_no_ref.nurbs = nrbdegelev(geom_no_ref.nurbs,[1 0]); % this way I get degree two in each direction
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
        
    case 3
        % Physical domain, defined as NURBS map given in a text file
        problem_data.geo_name = 'geo_thick_ring.txt';
        % initial import of geometry
        geom_no_ref  = geo_load (problem_data.geo_name);
% %         geom_no_ref.nurbs = nrbdegelev(geom_no_ref.nurbs,[1 0 1]); % this way I get degree two in each direction
% %         problem_data.D = geom_no_ref.rdim;
        
        % Type of boundary conditions for each side of the domain
        problem_data.nmnn_sides   = [];
        problem_data.drchlt_sides = [ 1 2 3 4 5 6];
        
        % Physical parameters
        problem_data.c_diff  = @(x, y, z) ones(size(x));
        
        % Source and boundary terms
        problem_data.f = @(x,y,z) (20*x.^3.*y.^2-30*x.*y.^2+12*x.*y.^4+2*x.^5+30*x.*y.^4-10*x.^3-60*x.*y.^2+24*x.^3.*y.^2+8*x).*sin(pi*z) - pi^2*(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.^2.*sin(pi*z);
        problem_data.g = @(x, y, z, ind) zeros(size(x));
        problem_data.h = @(x, y, z, ind) zeros(size(x));
        
        problem_data.uex = @(x,y,z) -(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.^2.*sin(pi*z);
        problem_data.dx_uex = @(x,y,z) -(4*x.^3-10*x+4*x.*y.^2).*x.*y.^2.*sin(pi*z) - (x.^4+y.^4-5*x.^2-5*y.^2+2*x.^2.*y.^2+4).*y.^2.*sin(pi*z);
        problem_data.dy_uex = @(x,y,z) -(4*y.^3-10*y+4*y.*x.^2).*x.*y.^2.*sin(pi*z) - (x.^4+y.^4-5*x.^2-5*y.^2+2*x.^2.*y.^2+4).*(2*y).*x.*sin(pi*z);
        problem_data.dz_uex = @(x,y,z) -(x.^2+y.^2-1).*(x.^2+y.^2-4).*x.*y.^2.*pi.*cos(pi*z);
        problem_data.graduex = @(x,y,z) cat (1,  reshape (problem_data.dx_uex(x,y,z), [1, size(x)]), ...
            reshape (problem_data.dy_uex(x,y,z), [1, size(x)]),...
            reshape (problem_data.dz_uex(x,y,z), [1, size(x)]));
        
end

if (isfield (geom_no_ref, 'nurbs'))
  ndim = numel (geom_no_ref.nurbs.order);
end

% % D = problem_data.D;                   
                   
% discretization parameters (p and h)

method_data.degree = 3*ones(1,ndim); % Degree of the splines, obtained by k-refinement of geometry, 
                                     % e.g [3 3] or [3 3 3]
method_data.nsub   = 8*ones(1,ndim); % will divide each subinterval of the original knot span in nsub many subinterval
                                     % e.g. [8 8] or [8 8 8]
                                     % i.e., we add nsub-1 knots in each interval of the original knotline.
                                     % note that if the original geometry has even number of subintervals, all possible
                                     % refinements with this strategy will have even subintervals

%% ================================  geometry ================================


% % % make sure degree is identical in each direction
% % if length(unique(geom_no_ref.nurbs.order))>1
% %     error('different degrees in different dim')
% % end

% Compute how many degrees to elev.
% % compute how many degrees to elev. If 0 or less, something is wrong
degelev  = max (method_data.degree - (geom_no_ref.nurbs.order-1), 0);
if (any (degelev < 0))
% %     error(strcat('max(degelev)=',num2str(max(degelev)),', stop! You need to fix the user set simulation degree'))
    warning('The degree provided is lower than the degree of the original geometry')
elseif (any (method_data.degree < 2))
    error ('Collocation for the Laplacian requires at least C^1 continuity. Degree should be at least 2')
end
nurbs_ptemp = nrbdegelev (geom_no_ref.nurbs, degelev);

% next h-ref (so the overall procedure is k-ref). 
[rknots, zeta, new_knots] = kntrefine (nurbs_ptemp.knots, method_data.nsub-1, nurbs_ptemp.order-1, nurbs_ptemp.order-2);

nurbs = nrbkntins (nurbs_ptemp, new_knots);
geo_refined = geo_load (nurbs);

% % for d = 1:ndim
% %         
% %     knotline = unique(nurbs_ptemp.knots{d});
% %     nb_sub0 = length(knotline)-1;
% %     NN = method_data.n_sub(d);
% %     new_knots{d} = [];
% % 
% %     for i = 1:nb_sub0
% %         tmp_new = linspace(knotline(i),knotline(i+1),NN+1);
% %         tmp_new(1)=[]; tmp_new(end)=[];
% %         new_knots{d}(end+1:end+method_data.n_sub(d)-1) = tmp_new;
% %         clear tmp_new
% %     end
% %     
% % end
% % 
% % % this is the final domain to use for computations
% % nurbs = nrbkntins( nurbs_ptemp, new_knots);

% these are the knot lines (with repetitions)
knots = geo_refined.nurbs.knots;


%% ================================ Collocation ================================

% Now we need to generate the collocation points. Since GeoPDEs is based on
%  a mesh structure, we generate an auxiliary "mesh" with only one point 
%  per element, using the midpoints between knots. 
%  Although this is not the most efficient way to implement
%  collocation, it allows us to use all the functionality in GeoPDEs, in
%  particular the evaluation of basis functions, without further changes.

% % % now we need to generate the collocation points. We define a sequence of coll_pts and create
% % % an auxiliary knot line so that we have 1 point per element. Then we build a mesh object
% % % over this knot line and the collocation points. Then GeoPDEs treats them as quadrature points, 
% % % but provides also a mean of evaluating the basis functions in those points, we go for it


%pts_case = 1; % equispaced
pts_case = 2; % greville -----------------------------------------------------------------> with greville I have coll-pts = DoFs so in 
                                                                                        % principle I do not need to use least squares. 
                                                                                        % However, I keep using it for the moment because 
                                                                                        % I want to be general in the choice of points
                                                                                        % and also the Dir BC are treated eliminating 
                                                                                        % the DoFs (matrix columns) but not rows
% pts_case = 3; % any external function prescribing points
coll_pts = cell(1,ndim);
switch pts_case
    case 1
        for d=1:ndim
% XXXX WHY 40? This should probably be in the input data.
            tmp = linspace(0,1,40); tmp(1)=[]; tmp(end)=[];
            coll_pts{d} = tmp; 
            clear tmp
        end
    case 2
        %  aveknt(knot_line,k)  returns the averages of successive  k-1  knots. To comply
        % with the definition of greville as g_j = ( zeta_{j+1} + ... + zeta_{j+p})/p
        % for an open knot line with n+p+1 and degree p, we then pass as second argument p+1
        for d = 1:ndim
% %             coll_pts{d} = aveknt(knots{d},method_data.degree(d)+1); 
            coll_pts{d} = aveknt (knots{d}, nurbs.order(d)); 
        end
    case 3
        % put your function here
        % coll_pts{d} = ...
end

% % % For each dim, we generate an auxiliary knot line, so that there is one coll point per element.
% % % we define this by taking as knot line [0 midpoint(coll-pt1, coll-pt2) midpoint(coll-pt2, coll-pt3) ... 1]
brk = cell(1,ndim);
for d = 1:ndim
    coll_pts{d} = coll_pts{d}(:)';
    if (numel(coll_pts{d}) > 1)
        brk{d} = [knots{d}(1), coll_pts{d}(1:end-1) + diff(coll_pts{d})/2, knots{d}(end)];
    else
        brk{d} = [knots{d}(1) knots{d}(end)];
    end
end

% This is the auxiliary mesh with one coll pt per element
coll_msh = msh_cartesian (brk, coll_pts, [], geo_refined, 'boundary', true,'der2',true);

% Now we generate the spline basis function set
space = sp_nurbs (geo_refined.nurbs, coll_msh);


% Now we evaluate each spline in the collocation points and collect everything into the design matrix, i.e. a
% matrix whose n-th row is the equation collocated in the n-th point. Remember that we have built coll_msh
% to have 1 pt per element, so the number of coll pts is exactly the number of elements of coll_msh

tot_nb_coll_pts = coll_msh.nel;
nb_dofs = space.ndof;
% % A = zeros (tot_nb_coll_pts, nb_dofs);
A = sparse (tot_nb_coll_pts, nb_dofs);

% Evaluate parameterization and basis functions
coll_mesh_eval = msh_precompute (coll_msh);
sp_evals  = sp_precompute (space, coll_mesh_eval, 'gradient', false, 'laplacian', true);

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

% Generate rhs. Recover coordinates of collocation points and evaluate forcing 
x = cell (1,ndim);
for d = 1:ndim
    x{d} = reshape (coll_mesh_eval.geo_map(d,:,:), coll_mesh_eval.nqn*coll_msh.nel, 1);
end
rhs = problem_data.f(x{:});

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

u_coll = zeros (space.ndof,1);
u_coll(internal_dofs) = AA\rr;


%% ================================ Galerkin for comparison ================================

% Even if we want to use the same space as before, we need to rebuild the
% space object because now mesh is different
nquad            = nurbs.order;     % Points for the Gaussian quadrature rule, equal to the order (degree + 1)
[gal_qn, gal_qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (nquad));
gal_msh          = msh_cartesian (knots, gal_qn, gal_qw, geo_refined);
gal_space = space.constructor (gal_msh);
% % gal_space = sp_nurbs(geo_refined.nurbs, gal_msh);
% % % gal_space    = sp_bspline (knots, method_data.degree, gal_msh);

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


% plot of solution (only for D=2)
if ndim==2
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
end


% compute errors of coll and gal. Create yet another mesh, which provides quad points for the error
quad_err_nquad      = 9*ones(1,ndim);     % Points for the Gaussian quadrature rule, like [9 9] or [9 9 9]
[quad_err_qn, quad_err_qw] = msh_set_quad_nodes (knots, msh_gauss_nodes(quad_err_nquad));
quad_err_msh      = msh_cartesian (knots, quad_err_qn, quad_err_qw, geo_refined);

% I also build another space, to evaluate the basis functions in the right quadrature points for error
quad_err_space = sp_nurbs( geo_refined.nurbs, quad_err_msh ); 


% Display errors of the computed solution in the L2 and H1 norm
format longg
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
 
