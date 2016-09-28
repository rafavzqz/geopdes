function Stiff_matrix = Stiff_WQ(msh, space, geometry, c_diff)
% This function generates the stiffness matrix using WQ

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
sp = space.constructor (new_msh);

% compute the function values and gradients at the quadrature points
for idim = 1:msh.ndim
    sp1d = sp.sp_univ(idim);
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
aux_val = zeros([msh.ndim msh.ndim prod(aux_size)]);
jac = msh.map_der(qn);
for i = 1:prod(aux_size)
	aux_val(:,:,i) = coeff(i)*abs(det(jac(:,:,i)))*inv(jac(:,:,i)'*jac(:,:,i));
end
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

for i = 1:N_dof

	ind = indices(i,:);
	for l = 1:d
		points{l} = Quad_rules(l).ind_points{ind(l)}; 
		j_act{l} = Connectivity(l).neighbors{ind(l)}; 
		len_j_act(l) = length(j_act{l});
	end
	i_nonzeros = prod(len_j_act);
	
	for k1 = 1:d % derivative index on the test function
		for k2 = 1:d % derivative index on the trial function
			C = aux_val(points{:},k1,k2); % coefficient tensor
			for l = d:-1:1 % sum_factorization loop
				if k1 == k2 && l == k1
					Q = Quad_rules(l).quad_weights_11{ind(l)};
					B = BSder{l,ind(l)}(1:len_j_act(l),:);
				elseif k1 ~= k2 && l == k1
					Q = Quad_rules(l).quad_weights_10{ind(l)};
					B = BSval{l,ind(l)}(1:len_j_act(l),:);				
				elseif k1 ~= k2 && l == k2
					Q = Quad_rules(l).quad_weights_01{ind(l)};
					B = BSder{l,ind(l)}(1:len_j_act(l),:);
				else
					Q = Quad_rules(l).quad_weights_00{ind(l)};
					B = BSval{l,ind(l)}(1:len_j_act(l),:);
				end
				B = bsxfun(@times,Q,B);
				if d == 3
					C = tprod(B,C,l);
				elseif d == 2
					if l == 1
						C = B*C;
					elseif l == 2
						C = C*B';
					end
				end
			end
			values(ncounter+1:ncounter+i_nonzeros) = values(ncounter+1:ncounter+i_nonzeros) + C(:)';
		end
	end
	
	% row indices
	rows(ncounter+1:ncounter+i_nonzeros) = i;
	
	% compute the column indices
	i_col = zeros(d,i_nonzeros);
	for l = 1:d
		rep = len_j_act; rep(l) = 1;
		perm = ones(1,d); perm(l) = len_j_act(l);
		ap = repmat(reshape(j_act{l}',perm),rep);
		i_col(l,:) = ap(:)';
	end
	cols(ncounter+1:ncounter+i_nonzeros) = 1 + (n_size.^(0:d-1))*(i_col-1);
% 	ap = cell(1,d);
% 	[ap{:}] = ndgrid(j_act{:});
% 	for k = d:-1:2
% 		cols(ncounter+1:ncounter+i_nonzeros) = cols(ncounter+1:ncounter+i_nonzeros) + (ap{k}(:)' - 1)*(n_size(k)^(k-1));
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

