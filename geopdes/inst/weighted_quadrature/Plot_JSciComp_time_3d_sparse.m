function [Time_NEW,Time_GeoPdes] = Plot_JSciComp_time_3d_sparse(pp1,nnel1)

if (nargin==0)
  pp1 = 2:8;
  nnel1 = 10;
end

problem_data.geo_name = 'geo_cube.txt';
% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];
% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y, z)  zeros (size (x));
problem_data.g = @test_cube_g_nmnn;
problem_data.h = @(x, y, z, ind) exp (x + z) .* sin (y);

for p = pp1
  q = p-1; % q � la regolarit� delle funzioni
  for nel = nnel1

    degree     = [p p p];       % Degree of the splines
    regularity = [p-1 p-1 p-1];       % Regularity of the splines
    nsub       = [nel nel nel];       % Number of subdivisions
    nquad      = [p+1 p+1 p+1];       % Points for the Gaussian quadrature rule
    geometry   = geo_load ('geo_cube.txt');
    [knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
    rule     = msh_gauss_nodes (nquad);
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh      = msh_cartesian (zeta, qn, qw, geometry);
    space    = sp_bspline (knots, degree, msh);

    disp(sprintf('Caso p=%0.5g, nel=%0.5g',p,nel))        
    dist_knots=linspace(0,1,nel+1); %zeta
    Knot_vect=augknt(dist_knots,p+1,p-q); %knots
    time_costruz_1D = clock;
    [my_JJ,p,dist_knots,mult_kv]=Index1D(Knot_vect);
    [my_QQ, Cond_G]=global_mid(p,q,dist_knots);
    
    sp1d = space.sp_univ(1);
%     for ii = 1:sp1d.ndof; my_JJ.neighbors2{ii} = unique (sp1d.connectivity(:,sp1d.supp{ii}))'; end
    my_JJ.neighbors = cellfun (@(x) unique (sp1d.connectivity(:,x)).', sp1d.supp, 'UniformOutput', false);
    my_JJ.num_neigh = cellfun (@numel, my_JJ.neighbors);
    my_JJ.first_elem = cellfun (@(x) x(1), sp1d.supp); % NOT USED
    my_JJ.last_elem = cellfun (@(x) x(end), sp1d.supp); % NOT USED
    
% Generate a new GeoPDEs mesh with the new quadrature points
    all_points = unique ([my_QQ.quad_points{:}]);
    my_QQ.all_points = all_points;
    for ii = 1:sp1d.ndof % Probably this can be improved
      [~,my_QQ.ind_points{ii}] = ismember (my_QQ.quad_points{ii}, all_points);
    end
    for jj = 1:3
      brk{jj} = [all_points(1), all_points(1) + diff(all_points)/2, all_points(end)];
      qn{jj} = all_points;
    end
    new_msh = msh_cartesian (brk, qn, [], geometry);
    sp = space.constructor (new_msh);
    sp1d = sp.sp_univ(1);

    % And we can probably improve also this one
    clear my_BBS
    for ii = 1:sp1d.ndof
      shape_funs = sp1d.shape_functions(:,:,my_QQ.ind_points{ii});
      conn = sp1d.connectivity(:,my_QQ.ind_points{ii});
      nsh = sp1d.nsh(my_QQ.ind_points{ii});
      for iq = 1:my_QQ.nquad_points(ii)
        [~,ind_supp,ind_elem] = intersect (my_JJ.neighbors{ii}, conn(:,iq));
        my_BBS{ii}(ind_supp,iq) = reshape (shape_funs(:,ind_elem,iq), [numel(ind_supp), 1]);
      end
    end
    
    BBS=BS_QP(p,q,dist_knots,Knot_vect,my_QQ,my_JJ);
    TEMPI_costruz_1D(p,nel)=etime(clock,time_costruz_1D);
    disp(sprintf('Costruzioni 1D %0.5g',TEMPI_costruz_1D(p,nel)))
    funz=@(x,y) 0.*kron(sin(x),ones(length(y),1))+0.*kron(ones(1,length(x)),y') +1;      
    total_time = tic;
    MM_3D_veloce = massSF_veloce3d(p,q,dist_knots,my_QQ,funz,my_JJ,BBS);
    TEMPI_MASS_3D_veloce(p,nel)=toc(total_time);
    disp(sprintf('Costruz matrice massa 3D Nuovo Metodo %0.5g',TEMPI_MASS_3D_veloce(p,nel)))        
    clear MM_3D_veloce

    end
 end

disp('Finito calcoli')
Time_NEW=TEMPI_MASS_3D_veloce(pp1,nnel1)+TEMPI_costruz_1D(pp1,nnel1);
disp('Tempi nuovo tutto include BS e costruz quadratura')
disp(Time_NEW)

% Time_GeoPdes=timeGeoPDES(pp1,nnel1)+TEMPI_Mesh_geopdes(pp1,nnel1);
% disp('Tempi GeoPdes tutto include geometria')
% disp(Time_GeoPdes)

end



function [my_JJ,p,dist_knots,mult_kv]=Index1D(Knot_vect)
% function [JJ,p,dist_knots,mult_kv]=Index1D_b(Knot_vect)
% Costruisce la matrice di connettivit` 1D
%
% JJ in output e' una matrice N_dof, (2p+1)+3
%       in JJ(ii,1) contiene il numero di funzioni con le quali la funzione phi(ii) si isnteseca
%       in JJ(ii,2:JJ(ii,1)+1) contiene gli indici delle funzioni con le quali la funzione phi(ii) si isnteseca 
%       in JJ(ii,end-1) contiene l'indice del primo elemento del supporto
%           della funzione phi(ii): n_1
%       in JJ(ii,end) contiene l'indice dell'ultimo elemento del supporto 
%           della funzione phi(ii): n_end
my_JJ = struct ('num_neigh', [], 'neighbors', [], 'first_elem', [], 'last_elem', []);

n_kv=max(size(Knot_vect));
dist_knots=Knot_vect(1);
kk=1;
mult_kv=0;

for jj=1:n_kv
    if Knot_vect(jj)>dist_knots(kk)
        dist_knots=[dist_knots,Knot_vect(jj)];
        mult_kv=[mult_kv,1];
        kk=kk+1;
    else
        mult_kv(kk)=mult_kv(kk)+1;
    end
end
p=mult_kv(1)-1;
n_dof=n_kv-(p+1);
n_el=length(dist_knots)-1;
% JJ=zeros(n_dof,(2*p+1)+3);
for ii=1:n_dof
    jjact=max([1,ii-p]):min([n_dof,ii+p]);
%     JJ(ii,1)=length(jjact); 
    my_JJ.num_neigh(ii) = numel (jjact);
%     JJ(ii,2:JJ(ii,1)+1)=jjact; 
    my_JJ.neighbors{ii} = jjact;
    for ll=1:n_el+1
       if  Knot_vect(ii) == dist_knots(ll)
%          JJ(ii,end-1)=ll; 
           my_JJ.first_elem(ii) = ll;
       end
       if  Knot_vect(ii+p+1) == dist_knots(ll)
%          JJ(ii,end)=ll-1; 
           my_JJ.last_elem(ii) = ll-1;
           break
       end
    end
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

function [my_QQ,Cond]=global_mid(p,q,dist_knots)

% QQ tensore di dimensione 3 x 2p+3 x N_dof che contiene informazione sui
% punti/pesi di quadratura
% QQ(3,1,ii) contiene il numero di punti di quadratura nel supporto della
% funzione di base ii 
% QQ(2,1:Loc_N_weights,ii) contiene i valori dei pesi di quadratura
% relativi alla funzione ii 
% QQ(1,1:Loc_N_weights,ii) contiene le ascisse dei punti di quadratura nel
% supporto della funzione di base ii 

n=length(dist_knots)-1;
N_dof=n*(p+1)-(n-1)*(q+1);
Knot_vect=augknt(dist_knots,p+1,p-q);
warning off
[xG_ref,wG_ref]=ggauss_Nab_xw(p+1);

[my_JJ,p,dist_knots_no,mult_kv]=Index1D(Knot_vect);
%%%% Prova per far funzionare spaziature qualsiasi
 %  Knot_vect(1)=Knot_vect(1)-100;
 %  Knot_vect(end)=Knot_vect(end)+100;
%%%% 
QQ=zeros(3,2*p+3,N_dof);
my_QQ = struct ('nquad_points', [], 'quad_weights', [], 'quad_points', []);

    n_el=length(dist_knots)-1;
    nodes=linspace(dist_knots(1),dist_knots(2),p+1);
    if n_el > 1
        if n_el==2
            nodes=nodes(1:end-1);
        end
        nodes=sort([nodes,linspace(dist_knots(end-1),dist_knots(end),p+1)]);
    end
    if n_el > 2
        nodes=sort([nodes,dist_knots(3:end-2)]);
        app=[];
        for iii=2:n_el-1
                app_1=linspace(dist_knots(iii),dist_knots(iii+1),3);
                app=[app,app_1(2:end-1)];
        end
        nodes=sort([nodes,app]);
    end

for ii=1:N_dof
    clear Loc*
    Loc_dist_knots(1)=Knot_vect(ii);
    kk=1;
    for jj=1:p+1
        if Knot_vect(ii+jj)>Loc_dist_knots(kk)
            Loc_dist_knots=[Loc_dist_knots;Knot_vect(ii+jj)];
            kk=kk+1;
        end
    end
    greville=[];
    for jjj=1:length(nodes)
        if Loc_dist_knots(1)<=nodes(jjj)
            if nodes(jjj)<=Loc_dist_knots(end)
                greville=[greville,nodes(jjj)];
            end
        end
     end
 
 %   Loc_Knot_vect=augknt(Loc_dist_knots,p+1,p-q);
    %%%%%%%%% trucco per far funzionare spcol %%%%%%%%%%
 %   Loc_Knot_vect(1)=Loc_Knot_vect(1)-100;
 %   Loc_Knot_vect(end)=Loc_Knot_vect(end)+100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %   Loc_xG=zeros(1,(p+1)*(kk-1)); Loc_wG=zeros(1,(p+1)*(kk-1));
    for ll=1:(kk-1)
        w=wG_ref.*(Loc_dist_knots(ll+1)-Loc_dist_knots(ll)); x=(Loc_dist_knots(ll+1)-Loc_dist_knots(ll)).*xG_ref + Loc_dist_knots(ll);
        Loc_xG((ll-1)*(p+1)+1:ll*(p+1))=x; Loc_wG((ll-1)*(p+1)+1:ll*(p+1))=w;
    end
    Loc_N_dof=(kk-1)*(p+1)-(kk-2)*(q+1);
    %QQ(3,1,ii)=Loc_N_dof;
    Loc_GAU=spcol(Knot_vect,p+1,Loc_xG);
%     Loc_GAU=Loc_GAU(:,JJ(ii,2:JJ(ii,1)+1));
    Loc_GAU=Loc_GAU(:,my_JJ.neighbors{ii});
    PHI_GG=spcol(Knot_vect(ii:ii+p+1),p+1,sort(Loc_xG));
    Loc_RHS=zeros(Loc_N_dof,1);
    for jj=1:Loc_N_dof
        Loc_RHS(jj)= sum(PHI_GG.*Loc_GAU(:,jj).*Loc_wG');
    end

    % Local Collocation matrix
    Loc_GG=spcol(Knot_vect,p+1, greville);    
    Loc_GG=Loc_GG(:,my_JJ.neighbors{ii});
        
    Cond(ii)=cond(Loc_GG);
    
    Loc_ww_greville=(Loc_GG'\Loc_RHS)';
    
    res = norm(Loc_GG'*Loc_ww_greville' - Loc_RHS)/norm(Loc_RHS);
    if res > 1e-12
        error('Too large residual')
    end
    
    Loc_N_weights=length(Loc_ww_greville);
    %QQ(1,1:Loc_N_dof,ii)=greville;
    QQ(1,1:Loc_N_weights,ii)=greville;
    my_QQ.quad_points{ii} = greville;
    %QQ(2,1:Loc_N_dof,ii)=Loc_ww_greville;
    QQ(2,1:Loc_N_weights,ii)=Loc_ww_greville;
    my_QQ.quad_weights{ii} = Loc_ww_greville;
    QQ(3,1,ii)=Loc_N_weights; %QQ(3,1,ii)=Loc_N_dof;
    my_QQ.nquad_points(ii) = Loc_N_weights;
end


end

function [BBS]=BS_QP(p,q,dist_knots,Knot_vect,my_QQ,my_JJ)

% BBS{ii} matrice che contiene i valori delle funzioni di
% base (nei punti di quadratura) che si intersecano con la funzione di base ii. 
% I punti di quadratura corrispondono alle colonne di BBS{ii}, mentre le funzioni
% corrispondono alle righe

n=length(dist_knots)-1;
N_dof=n*(p+1)-(n-1)*(q+1);
for ii=1:N_dof
%     jj_act=JJ1D(ii,2:JJ1D(ii,1)+1);
    jj_act = my_JJ.neighbors{ii};
%     Q=QQ(1,1:QQ(3,1,ii),ii);
    Q = my_QQ.quad_points{ii};
    for jjnd_act= 1:length(jj_act)%nq(m)
        BBS{ii}(jjnd_act,:)=spcol(Knot_vect(jj_act(jjnd_act):jj_act(jjnd_act)+p+1),p+1, Q)';
    end
end

end


function MM = massSF_veloce3d(p,q,dist_knots,my_QQ,funz,my_JJ,BBS)

% INPUT
% p           grado delle spline
% q           regolarit� delle spline
% dist_knots  vettore dei nodi senza ripetizioni
% QQ          tensore che contiene informazione sui punti/pesi di quadratura (vedi funzione global_mid)
% funz        input che contiene informazione sui coefficienti (attualmente inutilizzato)
% JJ1D        matrice che contiene informazione sulla connettivit� 1D (vedi funzione Index_1D)
% BBS         struttura che contiene i valori delle B-splines nei punti di quadratura (vedi funzione BS_QP)

d=3; % dimesnione
n=length(dist_knots)-1;
N_d=n*(p+1)-(n-1)*(q+1); % gradi di libert� per direzione
N_dof=(N_d)^d; % gradi di libert� totali

% nonzeros = sum(JJ1D(:,1))^3; % numero di nonzeri della matrice finale
nonzeros = sum (my_JJ.num_neigh.^3);
cols = zeros(1,nonzeros); rows = cols; values = cols;
ncounter = 0;

[ii1,ii2,ii3]=ind2sub([N_d,N_d,N_d], 1:N_dof);

for ii=1:N_dof

    i1 = ii1(ii); i2 = ii2(ii); i3 = ii3(ii);
% 	nq = QQ(3,1,[i1 i2 i3]); nq = nq(:)'; % nq vettore che contiene il numero di punti di quadratura relativi alla funzione [i1 i2 i3] in ogni direzione
	nq = my_QQ.nquad_points([i1 i2 i3]); nq = nq(:)'; % nq vettore che contiene il numero di punti di quadratura relativi alla funzione [i1 i2 i3] in ogni direzione
    
    %valutazione funz nei nodi!!!
    C = ones(nq);
    
    for ll=1:d

        if ll==1
%             Q = QQ(2,1:QQ(3,1,i3),i3);           % vettore dei pesi di quadratura relativi alla funzione i3
            Q = my_QQ.quad_weights{i3};           % vettore dei pesi di quadratura relativi alla funzione i3
%             jj_act_3 = JJ1D(i3,2:JJ1D(i3,1)+1);  % indici delle funzioni di base che si intersecano con la funzione i3
            jj_act_3 = my_JJ.neighbors{i3};      % indici delle funzioni di base che si intersecano con la funzione i3
            l_jj_act_3 = length(jj_act_3);
            BB = BBS{i3}(1:l_jj_act_3,:);        % matrice dei valori delle funzioni base (calcolate nei punti di quadratura) che si intersecano con la funzione i3
        elseif ll==2
%             Q = QQ(2,1:QQ(3,1,i2),i2);
            Q = my_QQ.quad_weights{i2};
%             jj_act_2 = JJ1D(i2,2:JJ1D(i2,1)+1);
            jj_act_2 = my_JJ.neighbors{i2};
            l_jj_act_2 = length(jj_act_2);
            BB = BBS{i2}(1:l_jj_act_2,:);
        elseif ll==3
%             Q = QQ(2,1:QQ(3,1,i1),i1);
            Q = my_QQ.quad_weights{i1};
%             jj_act_1 = JJ1D(i1,2:JJ1D(i1,1)+1);
            jj_act_1 = my_JJ.neighbors{i1};
            l_jj_act_1 = length(jj_act_1);
            BB = BBS{i1}(1:l_jj_act_1,:);
        end
        
        % Uso la funzione tprod per implementare il prodotto matrice-tensore
        % (riga 5 dello pseudocodice a pag. 14 di Calabr�,Sangalli,Tani)
        BB = bsxfun(@times,Q,BB);
        C = tprod(BB,C,d-ll+1);
        
    end

    C = C(:)';
    ii_nonzeros = length(C);
    
    app1 = repmat(jj_act_1',[1 l_jj_act_2 l_jj_act_3]);
    app2 = repmat(jj_act_2,[l_jj_act_1 1 l_jj_act_3]);
    app3 = repmat(reshape(jj_act_3,[1 1 l_jj_act_3]),[l_jj_act_1 l_jj_act_2 1]);
    ap1 = app1(:)'; ap2 = app2(:)'; ap3 = app3(:)';

    rows(ncounter+1:ncounter+ii_nonzeros) = ii;
    cols(ncounter+1:ncounter+ii_nonzeros) = (ap3-1)*N_d^2+(ap2-1)*N_d+ap1;
    values(ncounter+1:ncounter+ii_nonzeros) = C;
    ncounter = ncounter + ii_nonzeros;
   
end

MM = sparse (rows, cols, values, N_dof, N_dof); % assemblo la matrice


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
