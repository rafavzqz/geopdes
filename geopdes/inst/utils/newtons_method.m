% NEWTONS_METHOD: Solve a non-linear system of the form f(x) = 0. The method 
% stops when |f(x)| < tol or the maximum number of iterations is reached.
%
% USAGE:
%
%  [x,err,it] = newtons_method (f,x0,jac,tol,maxit)
%
% INPUT:
%
%  f: function handler of the non-linear system
%  
%  x0: starting vector
%
%  jac: function handler of the Jacobian of the non-linear system
%
%  tol: tolerance of the method
%
% OUTPUT:
%
%  x: solution such that f(x) = 0
%
%  err: residual of the last iteration
%
%  it: number of iterations
%
% Copyright (C) 2018 Luca Pegolotti
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [x,err,it] = newtons_method(f,x0,jac,tol,maxit)

res = f(x0);
err = norm(res);

it = 0;
x = x0;
while (err > tol && it < maxit)
    it = it + 1;
    disp(['Newton iteration ',num2str(it),' ...']);
    J = jac(x);
    a = J\res;
    x = x - a;
    res = f(x);
    err = norm(res);
    disp(['	done, error = ',num2str(err)]);
end



