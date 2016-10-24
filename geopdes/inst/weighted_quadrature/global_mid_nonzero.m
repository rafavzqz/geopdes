function Quad_rules = global_mid_nonzero (space)
% This function computes the 1D quadrature rules necessary to build the mass matrix 
% Note: still cannot handle non-uniform knot
%
% INPUT: 
% space: univariate spline space, as in space.sp_univ (see sp_scalar)
% 
% OUTPUT: structure Quad_rules, which contains the following fields
% Quad_rules.all_points    (array)               vector of quadrature points
% Quad_rules.nquad_points  (1 x ndof array)      number of quadrature points in the support of each basis function
% Quad_rules.quad_points   (1 x ndof cell-array) vector of quadrature points in the support of each basis function 
% Quad_rules.quad_weights  (1 x ndof cell-array) vector of quadrature weights in the support of each basis function 
% Quad_rules.ind_points    (1 x ndof cell-array) vector of indices of the points in quad_points{i} in the vector all_points 

knots = space.knots;
ndof = space.ndof;
degree = space.degree;
supp = space.supp;

Quad_rules = struct ('all_points',[],'nquad_points', [], 'quad_weights', [], 'quad_points', [],'ind_points',[]);
Quad_rules.nquad_points = zeros (ndof, 1);
Quad_rules.quad_weights = cell (ndof, 1);
Quad_rules.quad_points = cell (ndof, 1);
Quad_rules.ind_points = cell (ndof, 1);

% Quadrature points in the first, the last, and the internal elements
distinct_knots = unique (knots);
q_first_el = linspace (distinct_knots(1), distinct_knots(2), degree+2);
q_last_el = linspace (distinct_knots(end-1), distinct_knots(end), degree+2);
q_int_el = sort ([distinct_knots(3:end-2) 0.5*(distinct_knots(2:end-2)+distinct_knots(3:end-1))]);
all_points = [q_first_el q_int_el q_last_el];
all_points = unique (all_points);
Quad_rules.all_points = all_points;

% Construction of a global collocation basis, a matrix whose (i,j) entry 
%  is the value at the i-th quadrature point of the j-th basis function
iv = findspan (space.ndof-1, degree, all_points, knots);
B_loc = basisfun (iv, all_points, degree, knots);
num = numbasisfun (iv, all_points, degree, knots) + 1;
rows = repmat ((1:numel(all_points)).', 1, degree+1);
B_global = sparse (rows, num, B_loc, numel(all_points), ndof);

% Construction of a global collocation basis, with standard Gaussian quadrature points
%  It is assumed that the space is constructed on the standard Gaussian
%  quadrature. Otherwise, we have to evaluate the basis function on x_gauss, below.
nqn = size (space.shape_functions, 1);
nel = size (space.shape_functions, 3);
rows = repmat (reshape (1:nqn*nel, [nqn, 1, nel]), [1, space.nsh_max, 1]);
cols = repmat (reshape (space.connectivity, [1, space.nsh_max, nel]), [nqn, 1, 1]);
B_gauss = sparse (rows(:), cols(:), space.shape_functions(:), nqn*nel, ndof);

% Computation of the quadrature points for each function
%  The first (resp. last) function use p+1 uniform points in the first (last) element
%  Functions with support in the first (last) element also use these points.
% For internal functions, we remove the first and last points, where the function is zero
for ii = 1:ndof

    local_knots = knots(ii:ii+degree+1);
	distinct_local_knots = unique(local_knots);
	
    if (ii == 1)
        quad_points = all_points(1:degree+1);
    elseif (ii == ndof)
        quad_points = all_points(end-degree:end);
    else
        if (distinct_local_knots(1) == distinct_knots(1))
            quad_points = all_points(1:degree+2);
        else
            quad_points = [distinct_local_knots(1) 0.5*(distinct_local_knots(1)+distinct_local_knots(2)) distinct_local_knots(2)];
        end
        quad_points = [quad_points sort([distinct_local_knots(3:end-2) 0.5*(distinct_local_knots(2:end-2)+distinct_local_knots(3:end-1))])];

        if (distinct_local_knots(end) == distinct_knots(end))
            quad_points = [quad_points all_points(end-(degree+1):end)];
        else
            quad_points = [quad_points distinct_local_knots(end-1) 0.5*(distinct_local_knots(end-1)+distinct_local_knots(end)) distinct_local_knots(end)];
        end
        quad_points = quad_points(2:end-1); % Remove the first and lst knot, where the function vanishes
    end
    quad_points = unique (quad_points);
    Quad_rules.quad_points{ii} = quad_points;
    Quad_rules.nquad_points(ii) = numel (quad_points);
    [~,Quad_rules.ind_points{ii}] = ismember (quad_points, all_points);
end

gauss_rule = msh_gauss_nodes (degree + 1);
for ii = 1:ndof
    neighbors = unique (space.connectivity(:, supp{ii}));
    local_B = full (B_global(Quad_rules.ind_points{ii}, neighbors).');
	index_el_supp = supp{ii}'; 

    [qn, qw] = msh_set_quad_nodes (distinct_knots(index_el_supp(1):index_el_supp(end)+1), gauss_rule);
    x_gauss = qn(:).'; w_gauss = qw(:).';
    local_points = ((supp{ii}(1)-1)*nqn + 1) : supp{ii}(end)*nqn;
    local_rhs = B_gauss(local_points, neighbors).' * (B_gauss(local_points, ii).' .* w_gauss).';
    [Q,R] = qr(local_B',0);
    Quad_rules.quad_weights{ii} = (Q*(R'\local_rhs))';
end

end


% function [x,w]=ggauss_Nab_xw(N,a,b)
% % function [x,w]=ggauss_Nab_xw(N,a,b)
% %
% % Calculates the gauss nodes (x) and weights (w):
% % - N is the number of points of the quadrature rule
% % - [a,b] is the interval where the quadrature formula is requested.
% %
% % if only one input data is given, the formula is calculated in [0,1]
% %
% %
% % Uses the stable algorithm via closed form of recursive coefficients, for
% % details see FCACE2
% %
% 
% % DEF
% 
% ab=ones(N,2);% coef ricorsivi noti in forma chiusa
% ab(1,2)= 1;
% ab(2:N,2)= 1./(4*(4- 1./([1:N-1].^2) )); 
% ab(1:N,1)= ab(1:N,1)./2;
% J=zeros(N);
% for n=1:N, J(n,n)=ab(n,1); end
% for n=2:N
%   J(n,n-1)=sqrt(ab(n,2));
%   J(n-1,n)=J(n,n-1);
% end
% [V,D]=eig(J);
% [x,I]=sort(diag(D));
% V=V(:,I);
% w=ab(1,2)*V(1,:)'.^2;
% if nargin> 1
%  w=w.*(b-a); x=(b-a).*x + a;
% end
% 
% end
