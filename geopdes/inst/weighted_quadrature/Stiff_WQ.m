function Stiff_matrix = Stiff_WQ(msh, space, geometry, c_diff)
% This function generates the stiffness matrix using WQ

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

    Quad_rules(idim) = quadrule_stiff_fast(sp1d);
 
    brk{idim} = [space.knots{idim}(1), space.knots{idim}(end)];
    qn{idim} = Quad_rules(idim).all_points';
end

% Generate a new GeoPDEs mesh with the new quadrature points
new_msh = msh_cartesian (brk, qn, [], geometry);
space_wq = space.constructor (new_msh);

% compute the function values and gradients at the quadrature points
for idim = 1:msh.ndim
    sp1d = space_wq.sp_univ(idim);
    for ii = 1:sp1d.ndof
      BSval{idim,ii} = sp1d.shape_functions(Quad_rules(idim).ind_points{ii}, Connectivity(idim).neighbors{ii}).'; 
	  BSder{idim,ii} = sp1d.shape_function_gradients(Quad_rules(idim).ind_points{ii}, Connectivity(idim).neighbors{ii}).'; 
    end 	
end

% compute the coefficient at the quadrature points
x = cell (msh.ndim, 1);
[x{:}] = ndgrid(qn{:});
coeff = c_diff(x{:});
coeff = coeff(:);

aux_size = cellfun (@numel, qn);
jac = msh.map_der(qn);

C1 = reshape (coeff .* abs (geopdes_det__(jac)), 1, 1, new_msh.nqn, new_msh.nel);
jacT = permute (jac, [2 1 3 4]);

jacT_reshape = reshape (jacT, msh.ndim, msh.rdim, 1, new_msh.nqn, new_msh.nel);
jac_reshape = reshape (jac, 1, msh.rdim, msh.ndim, new_msh.nqn, new_msh.nel);
product = sum (bsxfun (@times, jacT_reshape, jac_reshape), 2);
product = reshape (product, msh.ndim, msh.ndim, new_msh.nqn, new_msh.nel);
aux_val = bsxfun (@times, C1, geopdes_inv__ (product));

aux_val = reshape(aux_val,[msh.ndim msh.ndim aux_size]);
aux_val = permute(aux_val,[3:2+msh.ndim 1 2]); % I put the first two indices in the last positions

Stiff_matrix = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,aux_val);

end


function Stiff_matrix = Stiff_fast(space,Quad_rules,Connectivity,BSval,BSder,aux_val)

d = ndims(aux_val)-2;
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

for ii = 1:N_dof

	ind = indices(ii,:);
	for ll = 1:d
		points{ll} = Quad_rules(ll).ind_points{ind(ll)}; 
		j_act{ll} = Connectivity(ll).neighbors{ind(ll)}; 
		len_j_act(ll) = length(j_act{ll});
	end
	i_nonzeros = prod(len_j_act);
	
	for k1 = 1:d % derivative index on the test function
		for k2 = 1:d % derivative index on the trial function
			C = aux_val(points{:},k1,k2); % coefficient tensor
			for ll = d:-1:1 % sum_factorization loop
				if (k1 == k2 && ll == k1)
					Q = Quad_rules(ll).quad_weights_11{ind(ll)};
					B = BSder{ll,ind(ll)}(1:len_j_act(ll),:);
				elseif (k1 ~= k2 && ll == k1)
					Q = Quad_rules(ll).quad_weights_10{ind(ll)};
					B = BSval{ll,ind(ll)}(1:len_j_act(ll),:);				
				elseif (k1 ~= k2 && ll == k2)
					Q = Quad_rules(ll).quad_weights_01{ind(ll)};
					B = BSder{ll,ind(ll)}(1:len_j_act(ll),:);
				else
					Q = Quad_rules(ll).quad_weights_00{ind(ll)};
					B = BSval{ll,ind(ll)}(1:len_j_act(ll),:);
				end
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
			values(ncounter+1:ncounter+i_nonzeros) = values(ncounter+1:ncounter+i_nonzeros) + C(:)';
		end
	end
	
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
% 	ap = cell(1,d);
% 	[ap{:}] = ndgrid(j_act{:});
% 	for k = d:-1:2
% 		cols(ncounter+1:ncounter+i_nonzeros) = cols(ncounter+1:ncounter+i_nonzeros) + (ap{k}(:)' - 1)*(n_index(k));
% 	end
% 	cols(ncounter+1:ncounter+i_nonzeros) = cols(ncounter+1:ncounter+i_nonzeros) + ap{1}(:)';
    
    ncounter = ncounter + i_nonzeros;
   
end

Stiff_matrix = sparse (rows, cols, values, N_dof, N_dof); % assemble the matrix

end


function Y = tprod(A,X,d)
% This function computes the product between matrix A and the 3D tensor X,
% in the direction d (d = 1,2 or 3)

[nx, ny, nz] = size(X);
m = size(A,1);

if d == 1
    
    Y = A*reshape(X,nx,ny*nz);
    Y = reshape(Y,m,ny,nz);
    
elseif d == 2
    
    Y = A*reshape(permute(X,[2 1 3]),ny,nx*nz);
    Y = permute(reshape(Y,m,nx,nz),[2 1 3]);
    
elseif d == 3
    
    Y = A*reshape(permute(X,[3 2 1]),nz,nx*ny);
    Y = permute(reshape(Y,m,ny,nx),[3 2 1]);
    
end

end
