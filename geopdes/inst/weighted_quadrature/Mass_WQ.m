function [Mass_matrix] = Mass_3D(msh, space, geometry, coeff_fun)

if (any (msh.nqn_dir ~= space.degree + 1))
  warning (['To compute the weights of the weighted quadrature the original ' ...
      'number of quadrature points should be p+1. Changing the points in the mesh.'])
  rule     = msh_gauss_nodes (space.degree+1);
  [qn, qw] = msh_set_quad_nodes (space.knots, rule);
  msh      = msh_cartesian (space.knots, qn, qw, geometry);
  space    = space.constructor (msh);
end

for idim = 1:msh.ndim
    sp1d = space.sp_univ(idim);
    Connectivity(idim).neighbors = cellfun (@(x) unique (sp1d.connectivity(:,x)).', sp1d.supp, 'UniformOutput', false);
    Connectivity(idim).num_neigh = cellfun (@numel, Connectivity(idim).neighbors);
    
%     Quad_rule(idim) = global_mid(sp1d);
    Quad_rule(idim) = global_mid_nonzero(sp1d);
 
    brk{idim} = [space.knots{idim}(1), space.knots{idim}(end)];
    qn{idim} = Quad_rule(idim).all_points';
end

% Generate a new GeoPDEs mesh with the new quadrature points
new_msh = msh_cartesian (brk, qn, [], geometry);
space_wq = space.constructor (new_msh);

for idim = 1:msh.ndim
    sp1d = space_wq.sp_univ(idim);
    for ii = 1:sp1d.ndof
      BSval{idim,ii} = sp1d.shape_functions(Quad_rule(idim).ind_points{ii}, Connectivity(idim).neighbors{ii}).';
    end 	
end

% calcolo le funzioni che compaiono nell'integrale
x = cell (msh.ndim, 1);
[x{:}] = ndgrid(qn{:});
coeff = coeff_fun(x{:});

aux_size = cellfun (@numel, qn);
coeff = reshape (coeff, aux_size);

jacdet = abs (geopdes_det__ (msh.map_der(qn)));
jacdet = reshape (jacdet, aux_size);
fun_val = jacdet.*coeff; 

Mass_matrix = Mass(space,Quad_rule,Connectivity,BSval,fun_val);

end


function Mass_matrix = Mass(space,Quad_rule,Connectivity,BSval,fun_val)

% INPUT
% space          spazio spline
% Quad_rule      tensore che contiene informazione sui punti/pesi di quadratura
% fun_val        input che contiene informazione sui coefficienti
% Connectivity   matrice che contiene informazione sulla connettivit� 1D
% BSval          struttura che contiene i valori delle B-splines nei punti di quadratura

d = ndims(fun_val);
N_dof = space.ndof;

nonzeros = prod(arrayfun(@(i)sum(Connectivity(i).num_neigh),1:d));
cols = zeros(1,nonzeros); rows = cols; values = cols;
ncounter = 0;

n_size = space.ndof_dir;
indices = cell(1,d);
[indices{:}] = ind2sub(n_size, 1:N_dof);
indices = cell2mat(indices);  indices = reshape(indices,[N_dof d]);
points = cell(1,d); j_act = cell(1,d); len_j_act = zeros(1,d); n_index = zeros(1,d);
for ll = 1:d
    n_index(ll) = prod(n_size(1:ll-1));
end

% row loop
for ii = 1:N_dof

	ind = indices(ii,:);
  
	for ll = 1:d
		points{ll} = Quad_rule(ll).ind_points{ind(ll)}; 
		j_act{ll} = Connectivity(ll).neighbors{ind(ll)}; 
		len_j_act(ll) = length(j_act{ll});
	end
	i_nonzeros = prod(len_j_act);
	C = fun_val(points{:}); % coefficient tensor
	for ll = d:-1:1 % sum_factorization loop
 		Q = Quad_rule(ll).quad_weights{ind(ll)};
 		B = BSval{ll,ind(ll)}(1:len_j_act(ll),:);
		B = bsxfun(@times,Q,B);
		if (d == 3)
			C = tprod(B,C,ll);
		elseif (d == 2)
			if (ll == 1)
				C = B*C;
			elseif (ll == 2)
				C = C*B';
			end
		end
	end
	values(ncounter+1:ncounter+i_nonzeros) = C(:)';

	% row indices
	rows(ncounter+1:ncounter+i_nonzeros) = ii;
	
	% compute the column indices
	i_col = zeros(d,i_nonzeros);
	for ll = 1:d
		rep = len_j_act; rep(ll) = 1;
		perm = ones(1,d); perm(ll) = len_j_act(ll);
		ap = repmat(reshape(j_act{ll}',perm),rep);
		i_col(ll,:) = ap(:)';      
	end
	cols(ncounter+1:ncounter+i_nonzeros) = 1 + n_index*(i_col-1);
  ncounter = ncounter + i_nonzeros;
   
end

Mass_matrix = sparse (rows, cols, values, N_dof, N_dof); % assemblo la matrice

end
