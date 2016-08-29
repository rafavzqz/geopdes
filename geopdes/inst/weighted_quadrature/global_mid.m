function Q = new_global_mid(space)
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

distinct_knots = unique(knots);
q_first_el = linspace(distinct_knots(1),distinct_knots(2),degree+1); % punti di quadratura nel primo elemento
q_last_el = linspace(distinct_knots(end-1),distinct_knots(end),degree+1);  % punti di quadratura nell'ultimo elemento
q_int_el = sort([distinct_knots(3:end-2) 0.5*(distinct_knots(2:end-2)+distinct_knots(3:end-1))]); % punti di quadratura negli elementi interni
all_points = [q_first_el q_int_el q_last_el];
Q.all_points = all_points;

for i = 1:ndof

    local_knots = knots(i:i+degree+1);
	distinct_local_knots = unique(local_knots);
	
    % calcolo i punti di quadratura per la funzione i
    if (i == 1) || (i == ndof) % gestisco a parte la prima e l'ultima funzione di base
        quad_points = linspace(distinct_local_knots(1),distinct_local_knots(2),degree+1);
    else
        if distinct_local_knots(1) == 0 % se la funzione ha supporto nel primo elemento, prendo i p+1 punti di quadratura corrispondenti
            quad_points = linspace(distinct_local_knots(1),distinct_local_knots(2),degree+1);
        else
            quad_points = [distinct_local_knots(1) 0.5*(distinct_local_knots(1)+distinct_local_knots(2))];
        end
        quad_points = [quad_points sort([distinct_local_knots(3:end-2) 0.5*(distinct_local_knots(2:end-2)+distinct_local_knots(3:end-1))])];
        if distinct_local_knots(end) == 1 % se la funzione ha supporto nell'ultimo elemento, prendo i p+1 punti di quadratura corrispondenti
            quad_points = [quad_points linspace(distinct_local_knots(end-1),distinct_local_knots(end),degree+1)];
        else
            quad_points = [quad_points distinct_local_knots(end-1) 0.5*(distinct_local_knots(end-1)+distinct_local_knots(end)) distinct_local_knots(end)];
        end
    end
    quad_points = unique(quad_points);
	Q.quad_points{i} = quad_points;
	Q.nquad_points(i) = length(quad_points);
	[~,Q.ind_points{i}] = ismember (quad_points, all_points);
	
	% costruisco la matrice di collocazione locale
    extended_local_knots = knots(max(1,i-degree):min(length(knots),i+2*degree+1)); % devo estendere il vettore dei nodi per beccare tutte le funzioni che mi servono
    local_B = spcol(extended_local_knots, degree+1, quad_points)';
    
	% costruisco il rhs con la quadratura gaussiana
	[x_guass_rif,w_gauss_rif]=ggauss_Nab_xw(degree+1,0,1);
    x_gauss=[];w_gauss=[];
	index_el_supp = supp{i}'; 
    for k = index_el_supp
        a =distinct_knots(k); b=distinct_knots(k+1);
        x_gauss = [x_gauss, x_guass_rif'*(b-a)+a];
        w_gauss = [w_gauss, w_gauss_rif'*(b-a)];
    end
    local_rhs = spcol(extended_local_knots,degree+1,x_gauss)'*((spcol(local_knots,degree+1,x_gauss)'.*w_gauss))';
	
	% calcolo i pesi
	quad_weights = local_B\local_rhs;
    Q.quad_weights{i} = quad_weights';
	
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

