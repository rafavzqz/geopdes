% mpi_gmres_LP: solve Ax=b with a left preconditioned parallel gmres
%                
%
% [sol, err] = mpi_gmres_LP (A, b, x0, LP, max_iter, comm)
%
% INPUT:
% 
% A:  matrix
% b:  right hand side
% x0: initial guess of the solution, can be used for restarting
% LP: left preconditioning matrix
% max_iter: maximum number of iterations, the memory used grows as
%           max_iter*n where n is the number of unknowns.
% comm: the mpi communication group
%
% OUTPUT:
%
% sol: approximative solution of the system
% err: residual that means err = b-A*sol
%
% Copyright (C) 2010 Andrea Bressan
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [sol, err] = mpi_gmres_LP (A, b, x0, LP, max_iter, comm)

n=length(b);

H=zeros(max_iter+1,max_iter);
V=zeros(n, max_iter+1);

v      = A*x0;
v      = MPI_Allreduce ( v, 'mpi_sum', comm);
v      = b - v;
v      = LP \ v;
v_norm = norm(v);
v      = v / v_norm;

V(:,1) = v;

for sol_iter=1:max_iter
  % compute the next basis element of the Krylov space w
  v = LP \ (A*v);
  v = MPI_Allreduce ( v, 'mpi_sum', comm);

  % ortogonalize w in respect to the previous base elements
  coef = V'*v;
  v    = v - V * coef;
  h    = norm(v,2);
  H(:,sol_iter) = coef';
  H(sol_iter+1,sol_iter) = h;
  if h==0
    break;
  end
  v = v / h;
  V(:,sol_iter+1) = v; 
end

HR = H(1:sol_iter+1, 1:sol_iter);
VR = V(:, 1:sol_iter);

aux = zeros(sol_iter+1,1);
aux(1,1) = v_norm;
y = zeros(sol_iter,1);
y = HR \ aux;
sol = x0 + VR*y;
approx = A*sol;
approx = MPI_Allreduce(approx, 'mpi_sum', comm);
err = b - approx;
end

%!test
%! dim=30
%! if ~MPI_Initialized()
%!    MPI_Init();
%! end
%! tag    = 33;
%! comm   = MPI_Comm_Load('ciao');
%! cpu_N  = MPI_Comm_size(comm);
%! cpu_ID = MPI_Comm_rank(comm);
%!
%! Am=rand(dim);
%! 
%! if cpu_ID==0
%!    b=rand(dim,1);
%!    x0=rand(dim,1);
%!    dest=setdiff(0:cpu_N-1,cpu_ID);
%!    MPI_Send(b,dest,tag,comm);
%!    MPI_Send(x0,dest,tag,comm);
%! else
%!    b=MPI_Recv(0,tag,comm);
%!    x0=MPI_Recv(0,tag,comm);
%! end
%! A=MPI_Allreduce(Am,'mpi_sum',comm);
%! [L,U]=lu(A);
%! [sol err]=mpi_gmres_LP(Am,b,x0, L,dim,comm);
%! sol_serial=A\b;
%! if(cpu_ID==0)
%!    disp(cat(2,sol,sol_serial))
%!    assert(max(sol-sol_serial)<= 1e-10);
%! end
%! MPI_Finalize();
