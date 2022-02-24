close all
clear all
warning ('off','geopdes:nrbmultipatch')

% PHYSICAL DATA OF THE PROBLEM
% 5 patches (L-shape)
% p1=[0 0]; p2=[1 0]; p3=[1 1]; p4=[-1/3 1/3]; p5=[0 -1]; p6=[-2/3 -1/3]; p7=[-1 -1]; p8=[-2/3 1/3]; p9=[-1 1];
% nrb(1) = nrb4surf(p1,p2,p4,p3);
% nrb(2) = nrb4surf(p1,p5,p6,p7);
% nrb(3) = nrb4surf(p6,p7,p8,p9);
% nrb(4) = nrb4surf(p1,p6,p4,p8);
% nrb(5) = nrb4surf(p8,p9,p4,p3);
% 
% problem_data.geo_name = nrb;

% % 3 patch surf
% coefs_1(:,3,:)=[0.125 -0.125 0.96875; 0.4625 0.1375 0.95; 0.8 0.4 0.2];
% coefs_1(:,2,:)=[-0.3125 0.1875 1.15625; 0.04375 0.44375 1.2625; 0.4 0.7 0.6];
% coefs_1(:,1,:)=[-0.75 0.5 0.1875;-0.375 0.75 0.5;  0.0 1.0 0.0];
% 
% coefs_2(:,3,:)=[0.125 -0.125 0.96875; 0.0125 -0.5125 0.9; -0.1 -0.9 0.18];
% coefs_2(:,2,:)=[0.4625 0.1375 0.95;  0.42375 -0.28125 1.134375; 0.385 -0.7 0.637];
% coefs_2(:,1,:)=[0.8 0.4 0.2; 0.835 -0.05 0.504;  0.87 -0.5 -0.0069];
% 
% coefs_3(:,3,:)=[0.125 -0.125 0.96875; -0.3125 0.1875 1.15625; -0.75 0.5 0.1875];
% coefs_3(:,2,:)=[0.0125 -0.5125 0.9; -0.39875 -0.25625 1.210625; -0.81 0.0 0.5975];
% coefs_3(:,1,:)=[-0.1 -0.9 0.18; -0.485 -0.7 0.463; -0.87 -0.5 -0.0069];
% 
% knots{1} = [0 0 0 1 1 1];
% knots{2} = [0 0 0 1 1 1];
% 
% butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
% butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
% butterf3 = nrbmak(permute(coefs_3,[3,2,1]),knots);
% 
% nrb(1)=butterf1;
% nrb(2)=butterf2;
% nrb(3)=butterf3;
% % nrb(1)=nrbtform (nrb(1), vecrot(pi/3,[1 2 3]));
% % nrb(2)=nrbtform (nrb(2), vecrot(pi/3,[1 2 3]));
% % nrb(3)=nrbtform (nrb(3), vecrot(pi/3,[1 2 3]));
% problem_data.geo_name = nrb;

%5 patch surf
coefs_1(:,1,:)=[0.125 -0.125 0.96875; 0.4375 0.2875 0.99375; 0.75 0.7 -0.0525];
coefs_1(:,2,:)=[-0.2375 0.3375 1.175; 0.06875 0.59375 1.0075; 0.375 0.85 0.3];
coefs_1(:,3,:)=[-0.6 0.8 0.0; -0.3 0.9 0.2; 0.0 1.0 0.0];

coefs_2(:,1,:)=[0.125 -0.125 0.96875; -0.2375 0.3375 1.175; -0.6 0.8 0.0];
coefs_2(:,2,:)=[-0.3875 -0.1625 1.0875; -0.58125 0.19625 0.88875; -0.775 0.555 0.182];
coefs_2(:,3,:)=[-0.9 -0.2 0.15; -0.925 0.055 0.207; -0.95 0.31 0.0014];

coefs_3(:,1,:)=[0.125 -0.125 0.96875; -0.3875 -0.1625 1.0875; -0.9 -0.2 0.15];
coefs_3(:,2,:)=[0.1125 -0.5458333333333333 0.8666666666666667; -0.29375 -0.5354166666666667 0.9264583333333334; -0.7 -0.525 0.38];
coefs_3(:,3,:)=[0.1 -0.9666666666666667 0.05555555555555555; -0.2 -0.9083333333333333 0.22833333333333333; -0.5 -0.85 0.0275];

coefs_4(:,1,:)=[0.125 -0.125 0.96875; 0.1125 -0.5458333333333333 0.8666666666666667; 0.1 -0.9666666666666667 0.05555555555555555];
coefs_4(:,2,:)=[0.5125 -0.2625 0.8375; 0.40625 -0.5854166666666667 0.6772916666666666; 0.3 -0.9083333333333333 0.12833333333333333];
coefs_4(:,3,:)=[0.9 -0.4 0.03; 0.7 -0.625 0.21; 0.5 -0.85 0.0275];

coefs_5(:,1,:)=[0.125 -0.125 0.96875; 0.5125 -0.2625 0.8375; 0.9 -0.4 0.03];
coefs_5(:,2,:)=[0.4375 0.2875 0.99375; 0.68125 0.12125 0.7625; 0.925 -0.045 0.269];
coefs_5(:,3,:)=[0.75 0.7 -0.0525; 0.85 0.505 0.0705; 0.95 0.31 0.0014];

knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 0 1 1 1];

butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
butterf3 = nrbmak(permute(coefs_3,[3,2,1]),knots);
butterf4 = nrbmak(permute(coefs_4,[3,2,1]),knots);
butterf5 = nrbmak(permute(coefs_5,[3,2,1]),knots);

nrb(1)=butterf1;
nrb(2)=butterf2;
nrb(3)=butterf3;
nrb(4)=butterf4;
nrb(5)=butterf5;
problem_data.geo_name = nrb;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [1 2 3 4 5 6 7 8 9 10];
problem_data.drchlt_sides = [];
problem_data.weak_drchlt_sides = [1 2 3 4 5 6 7 8 9 10];
        
% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));
% Source and boundary terms
C = 1;
p=-1;
% Bilaplacian
% problem_data.f = @(x,y,z) bilaplacian_rhs_Lshape2(x,y);
% % problem_data.g = @(x,y, ind) bilaplacian_Lshape_g_nmnn_r_8patches(x,y,ind);
% problem_data.g = @(x,y,z, ind) bilaplacian_Lshape_g_nmnn_r_5patches(x,y,ind); 
% problem_data.h = @(x, y, z, ind) solution_bilaplacian_Lshape (x, y);

% problem_data.f = @(x, y, z) (2./((1 + 4*x.^2 + 4*y.^2).^2)).*(cos(3*x).*((5 + 8*x.^4 + 38*y.^2 + 72*y.^4 + x.^2 .*(22 + 80*y.^2)).*cos(y) - 4*y.*(1 + 2*x.^2 + 2*y.^2).*sin(y)) + 12*x.*sin(3*x).*(-(1 + 2*x.^2 + 2*y.^2).*cos(y) + y.*(1 + 4*x.^2 + 4*y.^2).*sin(y)));
% problem_data.g = @(x, y, z, ind) (2*(3*x.*cos(y).*sin(3*x) + y.*cos(3*x).*sin(y)))./(1 + 4*x.^2 + 4*y.^2);
% problem_data.h = @(x, y, z, ind) cos(3*x) .* cos(y);
% problem_data.uex     = @(x, y, z) cos(3*x) .* cos(y);
% problem_data.graduex = @(x, y, z) cat (1, ...
%                        reshape ((-3*(1 + 4*y.^2).*cos(y).*sin(3*x) + 4*x.*y.*cos(3*x).*sin(y))./(1 + 4*x.^2 + 4*y.^2), [1, size(x)]), ...
%                        reshape ((12*x.*y.*cos(y).*sin(3*x) - (1 + 4*x.^2).*cos(3*x).*sin(y))./(1 + 4*x.^2 + 4*y.^2), [1, size(x)]), ...
%                        reshape ((2*(3*x.*cos(y).*sin(3*x) + y.*cos(3*x).*sin(y)))./(1 + 4*x.^2 + 4*y.^2), [1, size(x)]));  

problem_data.f = @(x, y, z) ( (4 * (1 + 2 * x.^2 + 2 * y.^2) ) ./ (1 + 4 * x.^2 + 4 * y.^2).^2);
problem_data.h = @(x, y, z, ind) z;
problem_data.g = @(x, y, z, ind) z;
problem_data.uex     = @(x, y, z) z;
problem_data.graduex = @(x, y, z) cat (1, ...
                       reshape (-( (2 * x) ./ (1 + 4 * x.^2 + 4 * y.^2) ), [1, size(x)]), ...
                       reshape (-( (2 * y) ./ (1 + 4 * x.^2 + 4 * y.^2) ), [1, size(x)]), ...
                       reshape ( (4 * (x.^2 + y.^2) ) ./ (1 + 4 * x.^2 + 4 * y.^2), [1, size(x)])); 

% Exact solution (optional)
problem_data.uex     = @(x, y, z) solution_bilaplacian_Lshape (x, y);
% problem_data.graduex = @(x, y, z) solution_bilaplacian_Lshape_grad (x, y);
% problem_data.hessuex = @(x, y, z) solution_bilaplacian_Lshape_hessian (x, y);

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];        % Degree of the splines
method_data.regularity  = [1 1];        % Regularity of the splines
method_data.nsub = [16 16];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);
method_data.space_type  = 'standard'; % 'simplified' (only children functions) or 'standard' (full basis)

% SOLVE AND DISPLAY RESULTS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u] =  mp_solve_bilaplace_C1 (problem_data, method_data);
if (isfield (problem_data, 'hessuex'))
    [err_h2, err_h1, err_l2, err_h2s, err_h1s] = sp_h2_error (hspace, hmsh, u, problem_data.uex, problem_data.graduex, problem_data.hessuex);
    if (plot_data.print_info) 
        fprintf('Error in H2 seminorm = %g\n', err_h2s); 
    end
end

% EXPORT VTK FILE
subplot(1,2,1)
npts = [51 51];
sp_plot_solution (u, hspace, geometry, npts); shading interp
vtk_pts = {linspace(0,1,npts(1)), linspace(0,1,npts(2))};
subplot(1,2,2)
for iptc = 1:hmsh.npatch
  F = reshape (geometry(iptc).map(vtk_pts), [hmsh.rdim, npts]);
  X = reshape (F(1,:), npts); Y = reshape (F(2,:), npts);
  surf(X, Y, problem_data.uex(X,Y));
  hold on; shading interp
end
warning ('on','geopdes:nrbmultipatch')
