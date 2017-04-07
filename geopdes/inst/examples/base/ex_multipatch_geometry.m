% EX_MULTIPATCH_GEOMETRY: generate and export to a file a multipatch geometry.

% The geometry is made of three patches: two rings, and a ruled surface
%  between an arc and a straight segment. The rings are built in two
%  different ways, to show some functions of the NURBS package

r1 = 0.2;
r2 = 0.5;
r3 = 1;
r4 = 2;

origin = [0 0 0];
start_angle = 0;
stop_angle = pi/6;

crv1 = nrbcirc (r1, origin, start_angle, stop_angle);
crv2 = nrbcirc (r2, origin, start_angle, stop_angle);
crv3 = nrbcirc (r3, origin, start_angle, stop_angle);

srf(1) = nrbruled (crv1, crv2);
srf(2) = nrbrevolve (nrbline([r2 0 0], [r3 0 0]), origin, [0 0 1], stop_angle);
srf(3) = nrbruled (crv3, nrbline([r4 0 0], [r4, r4*sin(stop_angle-start_angle)/cos(stop_angle-start_angle) 0]));

% Generate automatically the interface and boundary information
[interfaces, boundaries] = nrbmultipatch (srf);

% Using nrbmultipatch, one boundary is generated for each side of each patch
% To group several sides together, we must do it by hand
boundaries = [];

% Lower boundary
boundaries(1).nsides = 3;
boundaries(1).patches = [1 2 3];
boundaries(1).faces = [1 1 1];

% Upper boundary
boundaries(2).nsides = 3;
boundaries(2).patches = [1 2 3];
boundaries(2).faces = [2 2 2];

% Inner boundary
boundaries(3).nsides = 1;
boundaries(3).patches = 1;
boundaries(3).faces = 3;

% Outer boundary
boundaries(4).nsides = 1;
boundaries(4).patches = 3;
boundaries(4).faces = 4;

disp ('Exporting the geometry to the file geo_example_multipatch.txt')
nrbexport (srf, interfaces, boundaries, 'geo_example_multipatch.txt');

figure
for ii = 1:numel(srf); nrbkntplot(srf(ii)); hold on; end
