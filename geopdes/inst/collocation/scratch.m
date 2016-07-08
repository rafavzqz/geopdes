clear 

%% ================================ problem and choice of discretization  ================================

% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

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
problem_data.uex     = @(x, y) x(1-x) .* sin(pi*y);

                                     
% discretization parameters (p and h)

method_data.degree     = [3 3];       % Degree of the splines, obtained by k-refinement of geometry
method_data.h_level    = [3 3];       % 2^{hlev-1} nodes in each dir
            

%% ================================  geometry ================================



% initial import of geometry
geom_no_ref  = geo_load (problem_data.geo_name);
problem_data.D = geom_no_ref.rdim; 

% make sure degree is identical in each direction
degelev  = max (method_data.degree - (geom_no_ref.nurbs.order-1), 0);
if max(degelev)<=0,
    error(strcat('max(degelev)=',num2str(max(degelev)),', stop! You need to fix the user set simulation degree'))
end
% if so, to degree elev
nurbs_ptemp    = nrbdegelev (geom_no_ref.nurbs, degelev);

% next h-ref (so the overall procedure is k-ref). 
%  sidenote ------------------------------------------------------------------------> For the moment, exploiting diadic division is hard-coded for search of new points. See solve_laplace_iso for generic code
new_knots = cell(1,problem_data.D);
for d=1:problem_data.D
    new_knots{d}= ( 1 : (2^method_data.h_level(d)-1) ) /2^method_data.h_level(d) ;
end

% this is the final domain to use for computations
nurbs = nrbkntins( nurbs_ptemp, new_knots);
geo_refined = geo_load (nurbs);

% these are the knot lines (with repetitions)
knots=geo_refined.nurbs.knots;


%% ================================ Collocation ================================


% now we need to generate the collocation points. We generate a mesh object and place coll pts in the elements
% of the mesh. Geopdes treats them as quadrature points, but provides also a mean of evaluating the basis
% functions in those points, we go for it


% Let [x,w]=quadrule(n) be a quadrature formula over [-1 1]. Then I need its tensor version
% {[x1,w1],[x2,w2],[x3,w3]} = quadrulevect([n1 n2 n3])
% Consider
% knots = {knots1 knots2 knots3}
% where knotsK is the knot line along the K dir (even with repetition), defining N elements along K (total
% N1xN2xN3 elements).  Then
% [qn,qw] = msh_set_quad_nodes( knots, quadrulevect([n1 n2 n3]) )
% returns
% qn= {PT1 PT2 PT3}
% where PT1 is a matrix size n1 x N1 where each column are the points in each element


%[ qn, qw ] = msh_set_quad_nodes( knots, msh_gauss_nodes( [4 6] ) );
n_points_each_D = [3 3];
[ qn, qw ] = msh_set_quad_nodes( knots, n_points_quadrule( n_points_each_D ) );
msh = msh_cartesian( knots, qn, qw, geo_refined,'der2',true ); 


% now we generate the spline basis function set
space = sp_nurbs( geo_refined.nurbs, msh );

% now we evaluate each spline in the collocation points and collect everything into the design matrix, i.e. a
% matrix whose n-th row is the euqation collocated in the n-th point

tot_nb_coll_pts = prod(n_points_each_D)*msh.nel;  % sidenote ---------------------------------------------------> this one changes if alternating. How to force #pts = ndof?
nb_dofs = space.ndof;

A = zeros(tot_nb_coll_pts,nb_dofs);


% evaluate parameterization 
mesh_eval = msh_precompute(msh);
% evaluate basis fun 
sp_evals  = sp_precompute(space, mesh_eval, 'gradient', false, 'laplacian', true);

% sp_evals.shape_function_laplacian contains evaluations stored in a matrix of size
% [msh.nqn x nsh_max x msh.nel double]      
% i.e. for each node of an element (msh.nqn) only nonzero functions (nsh_max) 
% the nonzero functions are included in the matrix connectivity, whose size is (nsh_max x msh.nel vector)         
% so we need to unpack this information

for iel = 1:msh.nel
    
    %the nonzero funs on this element
    list_fun = sp_evals.connectivity(:,iel);
    
    % for each node, introduce a global numbering and store in A the evaluations
    for p = 1:mesh_eval.nqn        
        glob_p = p+ (iel-1)*mesh_eval.nqn;         
        A(glob_p,list_fun) = -sp_evals.shape_function_laplacians(p,:,iel); % sidenote ------------------------------> not this easy if instead we are solving -div{a grad u}=f
    end
    
end


% generate rhs. Recover coordinates of collocation points and evaluate forcing 
x={};
for idir = 1:mesh_eval.rdim
    x{idir} = reshape (mesh_eval.geo_map(idir,:,:), mesh_eval.nqn*msh.nel, 1);
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
%[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (gal_space, gal_msh, problem_data.h, problem_data.drchlt_sides); %-----------------------> this is not workgin???

gal_nb_boundaries = length(gal_space.boundary);
gal_boundary_dofs=[];
for bb=1:gal_nb_boundaries
    gal_boundary_dofs = union(gal_boundary_dofs,gal_space.boundary(bb).dofs);
end
gal_int_dofs = setdiff (1:gal_space.ndof, gal_boundary_dofs);

%gal_rhs(int_dofs) = gal_rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u_gal(gal_int_dofs) = gal_stiff_mat(gal_int_dofs, gal_int_dofs) \ gal_rhs(gal_int_dofs);




%% ============================== the comparison =======================================

plot(u_gal,'-x','DisplayName','Gal dofs')
hold on
plot(u_coll,'-or','DisplayName','Coll dofs')
legend show
