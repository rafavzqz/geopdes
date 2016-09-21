function Q = quadrule_stiff_fast(space)
% This function computes the 1D quadrature rules necessary to build the stiffness matrix 
% Note: still cannot handle non-uniform knot
%
% INPUT: 
% space           space of univariate splines
% 
% OUTUPT: structure Q
% Q.all_points        vector of quadrature points
% Q.nquad_points      vector of length length(all_points); nquad_points(i) = number of quadrature points in the support of the i-th basis function 
% Q.quad_points       structure array of length length(all_points); quad_points{i} = vector of quadrature points in the support of the i-th basis function 
% Q.ind_points        structure array of length length(all_points); ind_points{i} = vector of indices of the points in quad_points{i} in the vector all_points 
% Q.quad_weights_00   structure array of length length(all_points); quad_wights{i} = vector of quadrature weights relative to I00 in the support of the i-th function
% Q.quad_weights_10   structure array of length length(all_points); quad_wights{i} = vector of quadrature weights relative to I10 in the support of the i-th function
% Q.quad_weights_01   structure array of length length(all_points); quad_wights{i} = vector of quadrature weights relative to I01 in the support of the i-th function
% Q.quad_weights_11   structure array of length length(all_points); quad_wights{i} = vector of quadrature weights relative to I11 in the support of the i-th function

Q = struct ('all_points',[],'nquad_points', [], 'quad_weights_00', [], 'quad_weights_01', [], 'quad_weights_10', [], 'quad_weights_11', [], 'quad_points', [],'ind_points',[]);

knots = space.knots;
ndof = space.ndof;
degree = space.degree;

% Quadrature points in the first, the last, and the internal elements
%  I use eps to ensure that the values are inside the domain, and to avoid possible truncation errors
distinct_knots = unique (knots);
q_first_el = linspace (distinct_knots(1)+eps, distinct_knots(2), degree+2);
q_last_el = linspace (distinct_knots(end-1), distinct_knots(end)-eps, degree+2);
q_int_el = sort ([distinct_knots(3:end-2) 0.5*(distinct_knots(2:end-2)+distinct_knots(3:end-1))]);
all_points = [q_first_el q_int_el q_last_el];
all_points = unique (all_points);
Q.all_points = all_points;

% Computation of the quadrature points for each function. We consider only
% points were the function is nonzero.
for ii = 1:ndof
    if (ii == 1)
        quad_points = all_points(1:degree+1);
    elseif (ii == ndof)
        quad_points = all_points(end-degree:end);
    else
	    local_knots = knots(ii:ii+degree+1);
	    distinct_local_knots = unique(local_knots);
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
        quad_points = quad_points(2:end-1); % tolgo il primo e l'ultimo nodo, in cui la funzione vale zero 
    end
    quad_points = unique (quad_points);
	Q.quad_points{ii} = quad_points;
	Q.nquad_points(ii) = numel (quad_points);
	[~,Q.ind_points{ii}] = ismember (quad_points, all_points);
end

% Construction of a global collocation basis, a matrix whose (i,j) entry 
% is the value at the i-th quadrature point of the j-th basis function
iv = findspan (space.ndof-1, degree, all_points, knots); %%%%%%%%%%%%%%%%% ci va il -1 in ndof???
B_loc = basisfun (iv, all_points, degree, knots);
num = numbasisfun (iv, all_points, degree, knots) + 1;
rows = repmat ((1:numel(all_points)).', 1, degree+1);
B_global = sparse (rows, num, B_loc, numel(all_points), ndof);

% Construction of a global collocation basis, with standard Gaussian quadrature points
nqn = size (space.shape_functions, 1);
nel = size (space.shape_functions, 3);
rows = repmat (reshape (1:nqn*nel, [nqn, 1, nel]), [1, space.nsh_max, 1]);
cols = repmat (reshape (space.connectivity, [1, space.nsh_max, nel]), [nqn, 1, 1]);
B_gauss = sparse (rows(:), cols(:), space.shape_functions(:), nqn*nel, ndof);

% Construction of the derivative space. Note that we still use p+1 points
% for the gaussian quadrature (in this way the gaussian points are the same
% as in the original space and we do not need to re-compute the matrix B_gauss).
degree_der = degree - 1;
knots_der = knots(2:end-1);
ndof_der = ndof - 1;
nqn_der = degree_der + 2;
rule_der     = msh_gauss_nodes (nqn_der);
qn_der = msh_set_quad_nodes (distinct_knots, rule_der);
space_der = sp_bspline_1d_param(knots_der,degree_der,qn_der);

% Construction of a global basis for the derivative space, with given quadrature points
iv_der = findspan (ndof_der-1, degree_der, all_points, knots_der); %%%%%%%%%%%%%%%%% ci va il -1 in ndof???
B_loc_der = basisfun (iv_der, all_points, degree_der, knots_der);
num = numbasisfun (iv_der, all_points, degree_der, knots_der) + 1;
rows = repmat ((1:numel(all_points)).', 1, degree_der+1);
B_global_der = sparse (rows, num, B_loc_der, numel(all_points), ndof_der);

% Construction of a global basis for the derivative space, with standard Gaussian quadrature points
rows = repmat (reshape (1:nqn_der*nel, [nqn_der, 1, nel]), [1, space_der.nsh_max, 1]);
cols = repmat (reshape (space_der.connectivity, [1, space_der.nsh_max, nel]), [nqn_der, 1, 1]);
B_gauss_der = sparse (rows(:), cols(:), space_der.shape_functions(:), nqn_der*nel, ndof_der);

% Computation of the first derivatives of the basis functions at the Gaussian quadrature points
c = degree./(knots(degree+2:ndof+degree)-knots(2:ndof))';
B_prime_gauss = B_gauss_der*spdiags([-c c] ,[0 1],ndof-1,ndof);

Q.quad_weights_00 = quadrule_gen(space,space,Q.ind_points,B_global,B_gauss,B_gauss,degree+1);
Q.quad_weights_01 = quadrule_gen(space,space_der,Q.ind_points,B_global_der,B_gauss,B_gauss_der,degree+1);
Q.quad_weights_10 = quadrule_gen(space,space,Q.ind_points,B_global,B_prime_gauss,B_gauss,degree+1);
Q.quad_weights_11 = quadrule_gen(space,space_der,Q.ind_points,B_global_der,B_prime_gauss,B_gauss_der,degree+1);

end


function quad_weights = quadrule_gen(space_test,space_trial,ind_points,B_global_trial,B_gauss_test,B_gauss_trial,nqn_gauss)
% This function computes the weigthed quadrature rule for the test and trial function spaces in input, where the test functions are incorporated into the weights

distinct_knots = unique(space_test.knots);
ndof = space_test.ndof;
supp_test = space_test.supp;
connectivity_trial = space_trial.connectivity;

[x_gauss_rif,w_gauss_rif] = ggauss_Nab_xw(nqn_gauss,0,1); 

for i = 1:ndof
	% Construction of the local collocation matrix
    neighbors = unique (connectivity_trial(:, supp_test{i}));
    local_B = full (B_global_trial(ind_points{i}, neighbors).');
	% Construction of the rhs with gaussian quadrature
    x_gauss=[];w_gauss=[];
	index_el_supp = supp_test{i}'; 
    for k = index_el_supp
        a = distinct_knots(k); b = distinct_knots(k+1);
        x_gauss = [x_gauss, x_gauss_rif'*(b-a)+a];
        w_gauss = [w_gauss, w_gauss_rif'*(b-a)];
    end
    local_points = ((supp_test{i}(1)-1)*nqn_gauss + 1) : supp_test{i}(end)*nqn_gauss;
    local_rhs = B_gauss_trial(local_points, neighbors).' * (B_gauss_test(local_points, i).' .* w_gauss).';
    
	% computation of the weights
	quad_weights{i} = (local_B\local_rhs)';
end

end


function [x,w]=ggauss_Nab_xw(N,a,b)
% function [x,w]=ggauss_Nab_xw(N,a,b)
%
% Calculates the gauss nodes (x) and weights (w):
% - N is the number of points of the quadrature rule
% - [a,b] is the interval where the quadrature formula is requested.
%
% if only one input data is given, the formula is calculated in [0,1]
%
%
% Uses the stable algorithm via closed form of recursive coefficients, for
% details see FCACE2
%

% DEF

ab=ones(N,2);% coef ricorsivi noti in forma chiusa
ab(1,2)= 1;
ab(2:N,2)= 1./(4*(4- 1./([1:N-1].^2) )); 
ab(1:N,1)= ab(1:N,1)./2;
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[x,I]=sort(diag(D));
V=V(:,I);
w=ab(1,2)*V(1,:)'.^2;
if nargin> 1
 w=w.*(b-a); x=(b-a).*x + a;
end

end