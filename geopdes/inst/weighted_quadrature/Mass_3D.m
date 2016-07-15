function [Mass_matrix] = Mass_3D(msh, space, geometry, coeff_fun)

d = 3; 

for i = 1:d
    sp1d = space.sp_univ(i);
    Connectivity(i).neighbors = cellfun (@(x) unique (sp1d.connectivity(:,x)).', sp1d.supp, 'UniformOutput', false);
    Connectivity(i).num_neigh = cellfun (@numel, Connectivity(i).neighbors);
    
    Quad_rule(i) = global_mid(sp1d);
    %Quad_rule(i) = global_mid_nonzero(sp1d);
 
    brk{i} = [space.knots{i}(1), space.knots{i}(end)];
    qn{i} = Quad_rule(i).all_points';
end

% Generate a new GeoPDEs mesh with the new quadrature points
new_msh = msh_cartesian (brk, qn, [], geometry);
sp = space.constructor (new_msh);

for i = 1:d
    sp1d = sp.sp_univ(i);
    for ii = 1:sp1d.ndof
      BSval{i,ii} = sp1d.shape_functions(Quad_rule(i).ind_points{ii}, Connectivity(i).neighbors{ii}).';
    end 	
end

% calcolo le funzioni che compaiono nell'integrale
[x,y,z] = ndgrid(qn{1},qn{2},qn{3});
coeff = coeff_fun(x(:),y(:),z(:));
coeff = reshape (coeff, [length(qn{1}), length(qn{2}), length(qn{3})]);
jacdet = abs (geopdes_det__ (msh.map_der(qn)));
jacdet = reshape (jacdet, [length(qn{1}), length(qn{2}), length(qn{3})]);
fun_val = jacdet.*coeff; 

Mass_matrix = massSF_veloce3d(space,Quad_rule,Connectivity,BSval,fun_val);

end


function Mass_matrix = massSF_veloce3d(space,Quad_rule,Connectivity,BSval,fun_val)

% INPUT
% space          spazio spline
% Quad_rule      tensore che contiene informazione sui punti/pesi di quadratura
% fun_val        input che contiene informazione sui coefficienti
% Connectivity   matrice che contiene informazione sulla connettività 1D
% BSval          struttura che contiene i valori delle B-splines nei punti di quadratura

d=3; % dimesnione
N_dof = space.ndof;

% nonzeros = sum(JJ1D(:,1))^3; % numero di nonzeri della matrice finale
nonzeros = sum (Connectivity(1).num_neigh)*sum(Connectivity(2).num_neigh)*sum(Connectivity(3).num_neigh);
cols = zeros(1,nonzeros); rows = cols; values = cols;
ncounter = 0;

n1 = space.ndof_dir(1); n2 = space.ndof_dir(2); n3 = space.ndof_dir(3); 
[i1_vec,i2_vec,i3_vec]=ind2sub([n1,n2,n3], 1:N_dof);

for i=1:N_dof

    i1 = i1_vec(i); i2 = i2_vec(i); i3 = i3_vec(i);

    %valutazione funz nei nodi!
    C = fun_val(Quad_rule(1).ind_points{i1},Quad_rule(2).ind_points{i2},Quad_rule(3).ind_points{i3});
    
    for l=1:d

        if l == 1
            Q = Quad_rule(3).quad_weights{i3};           % vettore dei pesi di quadratura relativi alla funzione i3
            j_act_3 = Connectivity(3).neighbors{i3};      % indici delle funzioni di base che si intersecano con la funzione i3
            l_j_act_3 = length(j_act_3);
            B = BSval{3,i3}(1:l_j_act_3,:);        % matrice dei valori delle funzioni base (calcolate nei punti di quadratura) che si intersecano con la funzione i3
        elseif l == 2
            Q = Quad_rule(2).quad_weights{i2};
            j_act_2 = Connectivity(2).neighbors{i2};
            l_j_act_2 = length(j_act_2);
            B = BSval{2,i2}(1:l_j_act_2,:);
        elseif l == 3
            Q = Quad_rule(1).quad_weights{i1};
            j_act_1 = Connectivity(1).neighbors{i1};
            l_j_act_1 = length(j_act_1);
            B = BSval{1,i1}(1:l_j_act_1,:);
        end
        
        % Uso la funzione tprod per implementare il prodotto matrice-tensore
        B = bsxfun(@times,Q,B);
        C = tprod(B,C,d-l+1);
        
    end

    C = C(:)';
    i_nonzeros = length(C);
    
    app1 = repmat(j_act_1',[1 l_j_act_2 l_j_act_3]);
    app2 = repmat(j_act_2,[l_j_act_1 1 l_j_act_3]);
    app3 = repmat(reshape(j_act_3,[1 1 l_j_act_3]),[l_j_act_1 l_j_act_2 1]);
    ap1 = app1(:)'; ap2 = app2(:)'; ap3 = app3(:)';

    rows(ncounter+1:ncounter+i_nonzeros) = i;
    cols(ncounter+1:ncounter+i_nonzeros) = (ap3-1)*n3^2+(ap2-1)*n2+ap1;
    values(ncounter+1:ncounter+i_nonzeros) = C;
    ncounter = ncounter + i_nonzeros;
   
end

Mass_matrix = sparse (rows, cols, values, N_dof, N_dof); % assemblo la matrice


end


function Y = tprod(A,X,d)

% Questa funzione effettua il prodotto fra la matrice A e il tensore X
% nella direzione d = 1,2,3

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

