function Q = global_mid(space)
% Questa funzione genera le informazioni sulla quadratura 1D 
% N. B. Ancora non pu� gestire nodi ripetuti
%
% INPUT: 
% space           spazio di spline univariate
% 
% OUTUPT: struttura Q
% Q.all_points      vettore dei punti di quadratura
% Q.nquad_points    vettore di lunghezza ndof; nquad_points(i) = numero di punti di quadratura nel supporto della funzione i
% Q.quad_points     struttura di dimensione ndof; quad_points{i} = vettore dei punti di quadratura nel supporto della funzione i
% Q.quad_weights    struttura di dimensione ndof; quad_wights{i} = vettore dei punti di quadratura nel supporto della funzione i
% Q.ind_points      struttura di dimensione ndof; ind_points{i} = vettore degli indici dei dei punti di quad_points{i} nel vettore all_points

knots = space.knots;
ndof = space.ndof;
degree = space.degree;
supp = space.supp;

Q = struct ('all_points',[],'nquad_points', [], 'quad_weights', [], 'quad_points', [],'ind_points',[]);

% Quadrature points in the first, the last, and the internal elements
%  I use eps to ensure that the values are inside the domain, and to avoid possible truncation errors
distinct_knots = unique(knots);
q_first_el = linspace(distinct_knots(1)+eps,distinct_knots(2),degree+1); % punti di quadratura nel primo elemento
q_last_el = linspace(distinct_knots(end-1),distinct_knots(end)-eps,degree+1);  % punti di quadratura nell'ultimo elemento
q_int_el = sort([distinct_knots(3:end-2) 0.5*(distinct_knots(2:end-2)+distinct_knots(3:end-1))]); % punti di quadratura negli elementi interni
all_points = [q_first_el q_int_el q_last_el];
Q.all_points = all_points;

% Construction of a global collocation basis, with the value of the basis functions at the chosen points
iv = findspan (space.ndof, degree, all_points, knots);
B_loc = basisfun (iv, all_points, degree, knots);
num = numbasisfun (iv, all_points, degree, knots) + 1;
% B_global = sparse (numel(all_points), ndof);
% for ii = 1:numel(all_points)
%   B_global(ii, num(ii,:)) = B_loc(ii,:);
% end
rows = repmat ((1:numel(all_points)).', 1, degree+1);
B_global = sparse (rows, num, B_loc, numel(all_points), ndof);

% Construction of a global collocation basis, with standard Gaussian quadrature points
nqn = size (space.shape_functions, 1);
nel = size (space.shape_functions, 3);
% B_gauss = sparse (nqn * nel, ndof);
% for iel = 1:nel
%   B_gauss(nqn*(iel-1)+(1:nqn), space.connectivity(:,iel)) = space.shape_functions(:,:,iel);
% end
rows = repmat (reshape (1:nqn*nel, [nqn, 1, nel]), [1, space.nsh_max, 1]);
cols = repmat (reshape (space.connectivity, [1, space.nsh_max, nel]), [nqn, 1, 1]);
B_gauss = sparse (rows(:), cols(:), space.shape_functions(:), nqn*nel, ndof);

% Computation of the quadrature points for each function
%  The first (resp. last) function use p+1 uniform points in the first (last) element
%  Functions with support in the first (last) element also use these points.
for ii = 1:ndof

    local_knots = knots(ii:ii+degree+1);
	distinct_local_knots = unique(local_knots);
	
    % calcolo i punti di quadratura per la funzione i
    if (ii == 1)  % gestisco a parte la prima e l'ultima funzione di base
%         quad_points = linspace(distinct_local_knots(1),distinct_local_knots(2),degree+1);
        quad_points = all_points(1:degree+1);
    elseif (ii == ndof)
        quad_points = all_points(end-degree:end);
    else
        if (distinct_local_knots(1) == distinct_knots(1)) % se la funzione ha supporto nel primo elemento, prendo i p+1 punti di quadratura corrispondenti
%             quad_points = linspace(distinct_local_knots(1),distinct_local_knots(2),degree+1);
            quad_points = all_points(1:degree+2);
        else
            quad_points = [distinct_local_knots(1) 0.5*(distinct_local_knots(1)+distinct_local_knots(2))];
        end
        quad_points = [quad_points sort([distinct_local_knots(3:end-2) 0.5*(distinct_local_knots(2:end-2)+distinct_local_knots(3:end-1))])];
        if (distinct_local_knots(end) == distinct_knots(end)) % se la funzione ha supporto nell'ultimo elemento, prendo i p+1 punti di quadratura corrispondenti
%             quad_points = [quad_points linspace(distinct_local_knots(end-1),distinct_local_knots(end),degree+1)];
            quad_points = [quad_points all_points(end-(degree+1):end)];
        else
            quad_points = [quad_points distinct_local_knots(end-1) 0.5*(distinct_local_knots(end-1)+distinct_local_knots(end)) distinct_local_knots(end)];
        end
    end
    quad_points = unique (quad_points);
	Q.quad_points{ii} = quad_points;
	Q.nquad_points(ii) = length(quad_points);
	[~,Q.ind_points{ii}] = ismember (quad_points, all_points);
end

[x_gauss_rif,w_gauss_rif] = ggauss_Nab_xw(degree+1,0,1);
for ii = 1:ndof
	% costruisco la matrice di collocazione locale
%     extended_local_knots = knots(max(1,ii-degree):min(length(knots),ii+2*degree+1)); % devo estendere il vettore dei nodi per beccare tutte le funzioni che mi servono
%     local_B = spcol(extended_local_knots, degree+1, Q.quad_points{ii})';
%     local_knots = knots(ii:ii+degree+1);

%     neighbors = unique (num(Q.ind_points{ii},:));
    neighbors = unique (space.connectivity (:, supp{ii}));
    local_B = full (B_global(Q.ind_points{ii}, neighbors).');
	% costruisco il rhs con la quadratura gaussiana
	[x_guass_rif,w_gauss_rif]=ggauss_Nab_xw(degree+1,0,1);
    x_gauss=[];w_gauss=[];
	index_el_supp = supp{ii}'; 
    for k = index_el_supp
        a =distinct_knots(k); b=distinct_knots(k+1);
        x_gauss = [x_gauss, x_guass_rif'*(b-a)+a];
        w_gauss = [w_gauss, w_gauss_rif'*(b-a)];
    end
%     local_rhs = spcol(extended_local_knots,degree+1,x_gauss)'*((spcol(local_knots,degree+1,x_gauss)'.*w_gauss))';
    local_points = ((supp{ii}(1)-1)*nqn + 1) : supp{ii}(end)*nqn;
    local_rhs = B_gauss(local_points, neighbors).' * (B_gauss(local_points, ii).' .* w_gauss).';
	
	% calcolo i pesi
	quad_weights = local_B\local_rhs;
    Q.quad_weights{ii} = quad_weights';
	
end

end
function [x,w]=ggauss_Nab_xw(N,a,b)
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

