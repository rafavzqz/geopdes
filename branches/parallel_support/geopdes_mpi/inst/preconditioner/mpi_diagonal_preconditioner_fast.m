% MPI_DIAGONAL_PRECONDITIONER_FAST
%
% PL = mpi_diagonal_preconditioner_fast(A, Comm)
%
%
% INPUTS:
%
% A:    the matrix of the system memorized by the processor
% Comm: the MPI communication group 
%
% OUTPUTS:
%
% PL:   The left diagonal preconditioning. The coefficients are an estimate
%       of the norm of the rows of A. Each process i compute the sum
%       of the squared coefficients of each row of A_i  then the the result of
%       all processes are summed and the square root is calculated.
%       Computing the true norm requires summing the rows and then summing the
%       square of the coefficients, this is very expensive almost doubling the 
%       solving time because the complete rows of A must be sent.
%       This approach gives a slightly worse preconditioner for a very smaller
%       cost. Numerical test suggest that this is a good approach.
%
% Copyright (C) 2010 Andrea Bressan
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see  <http://www.gnu.org/licenses/>.

function PL = mpi_diagonal_preconditioner_fast (A, Comm)
  A2   = A.^2;
  RN   = sum ( A2, 2 );
  clear A2;
  D2   = MPI_Allreduce(RN, 'mpi_sum', Comm);
  PL   = sparse( diag( sqrt(D2) ) );
end