% mpi_Rgmres: solves Ax=b with a parallel restarted gmres
%
% [sol, err] = mpi_Rgmres (A, b, x0, rtol, max_iter,
%                         restart, comm, preconditioners, debug_print)
%
% INPUT:
%
% A : the matrix
% b : the right hand side
% x0: initial guess of the solution
% 
% rtol: gives the relative tollerance of the solution. The algorytm
%       stops whenever  the ratio between the norm of the residual
%       and that of b drops below rtol.
% max_iter: the max number of iterations
% restart: number of iterations between restarts
% comm:    the communication group
% preconditioners: a cell array of two matrices {PL,PR}: the left and
%                  right preconditioner respectively. Any of the two can
%                  be empty.
% debug_print: if 0 no info are printed, if 1 the iteration, residual norm
%              and the ratio between the b norm and the residual norm are
%              printed at every restart 
%
% OUTPUT:
%
% sol: the approximate solution
% err: the residual b-A*sol
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

function [sol, res] = mpi_Rgmres(A,b,x0, rtol, max_iter, restart, Comm,P, debug_print)

norm_res = Inf;
norm_b   = norm(b,2);
sol=x0;

args_num=nargin();

if args_num < 7 || args_num > 9
  error('mpi_Rgmres: requires between 7 to 9 arguments');
elseif args_num == 7
  P={[],[]};
  args_num=8;
elseif args_num == 9
  if isscalar(debug_print)
    if debug_print==0
      debug_print=@(bb, rres, nnormb, nnormres, iiter)0;
    else
      mpi_message(sprintf('Rgmres, restart: %d, max_iter: %d, rtol: %e, norm RHS: %e', restart, max_iter, rtol, norm_b),Comm);
      debug_print=@(bb, rres, nnormb,  nnormres, iiter)mpi_message(sprintf('Iter: %d;  ResNorm:  %e; RelNorm: %e;', iiter, nnormres, nnormres/nnormb ),Comm);
    end
  end
elseif args_num == 8
  debug_print=@(bb, rres, nnormb, nnormres,iiter)0;
end

if norm(sol,2)==0 % wrong initial guess of the solution
  if norm_b == 0  % ok the guess is the solution
    sol=zeros(size(x0));
    err=sol;
    return;
  else
    error('mpi_Rgmres: started with a null vector as starting point');
  end
end
iter=0;

if isempty(P{1}) && isempty(P{2})
				% NO PRECONDITIONING
  while (max_iter > iter && norm_res > rtol*norm_b)
    [sol,res]=mpi_gmres(A, b, sol, restart, Comm);
    norm_res=norm(res,2);
    iter=iter+restart;
    debug_print(b,res,norm_b,norm_res,iter);
  end
elseif isempty(P{1}) && ~isempty(P{2})
				% RIGHT PRECONDITIONING
  RP=P{2};
  while (max_iter > iter && norm_res > rtol*norm_b)
    [sol,res]=mpi_gmres_RP(A, b, sol, RP, restart, Comm);
    norm_res=norm(res,2);
    iter=iter+restart;
    debug_print(b,res,norm_b,norm_res,iter);
  end
elseif  ~isempty(P{1}) && isempty(P{2})
				% LEFT PRECONDITIONING
  LP=P{1};
  while (max_iter > iter && norm_res > rtol*norm_b)
    [sol,res]=mpi_gmres_LP(A, b, sol, LP, restart, Comm);
    norm_res=norm(res,2);
    iter=iter+restart;
    debug_print(b,res,norm_b,norm_res,iter);
  end
else
				% SPLIT PRECONDITIONING
  LP={P1};
  RP=P{2};
  while (max_iter > iter && norm_res > rtol*norm_b)
    [sol,res]=mpi_gmres_LRP(A, b, sol, LP, RP, restart, Comm);
    norm_res=norm(res,2);
    iter=iter+restart;
    debug_print(b,res,norm_b,norm_res,iter);
  end
end

end
