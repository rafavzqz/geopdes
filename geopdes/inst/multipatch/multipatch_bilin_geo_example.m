clear all

p1=[0 0]; p2=[1 0]; p3=[0 1]; p4=[1 1]; p5=[2,0]; p6=[2,1];
p7=[0 -1]; p8=[1 -1]; p9=[2 -1];
nrb(1) = nrb4surf(p2,p4,p1,p3);
nrb(2) = nrb4surf(p2,p5,p4,p6);
nrb(3) = nrb4surf(p2,p8,p5,p9);
nrb(4) = nrb4surf(p2,p1,p8,p7);
problem_data.geometry = nrb;

[geometry, boundaries, interfaces, subdomains, boundary_interfaces] = mp_geo_load (nrb);
[interfaces_all, vertices] = vertices_struct(boundaries, interfaces);