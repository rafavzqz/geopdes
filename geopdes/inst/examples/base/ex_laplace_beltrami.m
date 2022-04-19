% EX_LAPLACE_BELTRAMI: solve the Laplace-Beltrami equation on a NURBS surface, discretized with B-splines (non-isoparametric approach).

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS structure from the NURBS toolbox
% problem_data.geo_name = 'geo_roof.txt';

% % 4 patch surf hyperbolic smooth boundary
% coefs_1(:,1,:)=[-0.5 -0.5 0.0; -0.25 -0.5 -0.25; 0.0 -0.5 -0.25];
% coefs_1(:,2,:)=[-0.5 -0.2 0.3; -0.2375 -0.25 -0.0125; 0.025 -0.3 -0.05];
% coefs_1(:,3,:)=[-0.5 0.1 0.24; -0.225 0.0 -0.015; 0.05 -0.1 -0.0075];
% 
% coefs_2(:,1,:)=[0.5 -0.5 0.0; 0.5 -0.225 0.275; 0.5 0.05 0.2475];
% coefs_2(:,2,:)=[0.25 -0.5 -0.25; 0.2625 -0.2625 0.0; 0.275 -0.025 0.03];
% coefs_2(:,3,:)=[0.0 -0.5 -0.25; 0.025 -0.3 -0.05; 0.05 -0.1 -0.0075];
% 
% coefs_3(:,1,:)=[0.05 -0.1 -0.0075; 0.275 -0.025 0.03; 0.5 0.05 0.2475];
% coefs_3(:,2,:)=[0.0 0.2 0.0475; 0.25 0.2375 0.0125; 0.5 0.275 0.225];
% coefs_3(:,3,:)=[-0.05 0.5 -0.2475; 0.225 0.5 -0.275; 0.5 0.5 0.0];
% 
% coefs_4(:,1,:)=[-0.5 0.5 0.0; -0.5 0.3 0.2; -0.5 0.1 0.24];
% coefs_4(:,2,:)=[-0.275 0.5 -0.225; -0.25 0.25 0.0; -0.225 0.0 -0.015];
% coefs_4(:,3,:)=[-0.05 0.5 -0.2475; 0.0 0.2 0.0475; 0.05 -0.1 -0.0075];
% 
% knots{1} = [0 0 0 1 1 1];
% knots{2} = [0 0 0 1 1 1];
% 
% butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
% butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
% butterf3 = nrbmak(permute(coefs_3,[3,2,1]),knots);
% butterf4 = nrbmak(permute(coefs_4,[3,2,1]),knots);
% 
% nrb(1)=butterf1;
% % nrb(2)=butterf2;
% % nrb(3)=butterf3;
% % nrb(4)=butterf4;
% problem_data.geo_name = nrb;


% % 1 patch surf hyperbolic smooth boundary
coefs_1(:,1,:)=[-0.5 -0.5 0.0; 0.0 -0.5 -0.5; 0.5 -0.5 0.0];
coefs_1(:,2,:)=[-0.5 0.0 0.5; 0.0 0.0 0.0; 0.5 0.0 0.5];
coefs_1(:,3,:)=[-0.5 0.5 0.0; 0.0 0.5 -0.5; 0.5 0.5 0.0];

knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 0 1 1 1];

butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);

nrb(1)=butterf1;
problem_data.geo_name = nrb;


% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% % Source and boundary terms
% alp = 5; bet = 1; L = 1;
% g1 = @(x,y) (1-cos(atan(x./y))).*(1-sin(atan(x./y)));
% g2 = @(x,y) cos(atan(x./y)) + sin(atan(x./y)) - 4*cos(atan(x./y)).*sin(atan(x./y));
% gz = @(z) sin(pi*alp/L*z);
% dg1 = @(x,y) sin(atan(y./x)).*(1-sin(atan(y./x))) - cos(atan(y./x)).*(1-cos(atan(y./x)));
% 
% problem_data.f = @(x, y, z) bet*((alp*pi/L)^2 * g1(x,y) - g2(x,y)).*gz(z);
% problem_data.h = @(x, y, z, ind) bet * g1(x,y).*gz(z);
% 
% % Exact solution (optional)
% problem_data.uex     = @(x, y, z) bet * g1(x,y).*gz(z);
% problem_data.graduex = @(x, y, z) cat (1, ...
%                 reshape (-1./sqrt(x.^2+y.^2).*sin(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
%                 reshape ( 1./sqrt(x.^2+y.^2).*cos(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
%                 reshape ( alp*pi/L * cos(alp*pi/L*z) .* bet .* g1(x,y), [1, size(x)]));

% % z hyperbolic
% problem_data.f = @(x, y, z) z;
problem_data.f = @(x, y, z) ((8 .*(x.^2 - y.^2))./(1 + 4 *x.^2 + 4 *y.^2).^2);
% problem_data.g = @(x, y, z, ind) z;
problem_data.h = @(x, y, z, ind) z;
problem_data.uex     = @(x, y, z) z;
problem_data.graduex = @(x, y, z) cat (1, ...
                       reshape ((2*x)./(1 + 4*x.^2 + 4*y.^2), [1, size(x)]), ...
                       reshape (-((2*y)./(1 + 4*x.^2 + 4*y.^2)), [1, size(x)]), ...
                       reshape (((4*(x.^2 + y.^2))./(1 + 4*x.^2 + 4*y.^2)), [1, size(x)])); 
% problem_data.hessuex = @(x, y, z) zeros([3, 3, size(x)]);
% problem_data.lapuex  = @(x, y, z) -((8*(x.^2 - y.^2))./(1 + 4*x.^2 + 4*y.^2).^2);


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [1 1];       % Regularity of the splines
% method_data.nsub       = [9 9];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

h=[];

for i= 1:4
    fprintf ('Loop %d \n', i);
    method_data.nsub       = [2^(i+1) 2^(i+1)];      % Number of subdivisions
    % 3) CALL TO THE SOLVER
    
    % 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_laplace (problem_data, method_data);
    
    h = [h 1/sqrt(msh.nel(1))];
% Display errors of the computed solution in the L2 and H1 norm
[error_h1(i), error_l2(i)] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex)
end

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Laplace-Beltrami.vts';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,2,1)
surf (X, Y, eu)
title ('Numerical solution'), axis tight
subplot (1,2,2)
surf (X, Y, problem_data.uex (X,Y))
title ('Exact solution'), axis tight


figure
loglog(h, error_l2,'-o')
hold on
loglog(h, error_h1,'-o')
loglog(h, 100*h.^4,'-x')
loglog(h, 100*h.^3,'-x')
loglog(h, 10*h.^2,'-x')
legend("L^2 error", "H^1 error", "h^4", "h^3", "h^2",'Location','southeast');


% figure
% npts = [20 20];
% [eu, F] = sp_eval(u, space, geometry, vtk_pts, 'gradient');
% % F = reshape (geometry(iptc).map(vtk_pts), [hmsh.rdim, npts]);
% X = reshape (F(1,:), npts); Y = reshape (F(2,:), npts); Z = reshape (F(3,:), npts);
% gradex = problem_data.graduex(X,Y,Z);
% subplot (1,2,1)
% surf(X, Y, Z, squeeze(gradex(3,:,:)));
% shading interp
% subplot (1,2,2)
% surf(X, Y, Z, squeeze(eu(3,:,:)));
% shading interp


%!demo
%! ex_laplace_beltrami

%!test
%! problem_data.geo_name = 'geo_roof.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! alp = 5; bet = 1; L = 1;
%! g1 = @(x,y) (1-cos(atan(x./y))).*(1-sin(atan(x./y)));
%! g2 = @(x,y) cos(atan(x./y)) + sin(atan(x./y)) - 4*cos(atan(x./y)).*sin(atan(x./y));
%! gz = @(z) sin(pi*alp/L*z);
%! dg1 = @(x,y) sin(atan(y./x)).*(1-sin(atan(y./x))) - cos(atan(y./x)).*(1-cos(atan(y./x)));
%! problem_data.f = @(x, y, z) bet*((alp*pi/L)^2 * g1(x,y) - g2(x,y)).*gz(z);
%! problem_data.h = @(x, y, z, ind) bet * g1(x,y).*gz(z);
%! problem_data.uex     = @(x, y, z) bet * g1(x,y).*gz(z);
%! problem_data.graduex = @(x, y, z) cat (1, ...
%!                 reshape (-1./sqrt(x.^2+y.^2).*sin(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
%!                 reshape ( 1./sqrt(x.^2+y.^2).*cos(atan(y./x))*bet.*gz(z).*dg1(x,y), [1, size(x)]), ...
%!                 reshape ( alp*pi/L * cos(alp*pi/L*z) .* bet .* g1(x,y), [1, size(x)]));
%! method_data.degree     = [3 3];       % Degree of the splines
%! method_data.regularity = [2 2];       % Regularity of the splines
%! method_data.nsub       = [9 9];       % Number of subdivisions
%! method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_laplace (problem_data, method_data);
%! [error_h1, error_l2] = ...
%!            sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
%! assert (msh.nel, 81)
%! assert (space.ndof, 144)
%! assert (error_h1, 0.0407258650845198, 1e-15)
%! assert (error_l2, 0.00101798638120372, 1e-16)
