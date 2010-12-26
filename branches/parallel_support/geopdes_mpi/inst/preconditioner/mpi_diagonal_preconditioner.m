% MPI_DIAGONAL_PRECONDITIONER
%
% PL = mpi_diagonal_preconditioner(A, Comm)
%
%
% INPUTS:
%
% A:    the matrix of the system
% Comm: the MPI communication group 
%
% OUTPUTS:
%
% PL:   the diagonal preconditioning matrix whose entries are the
%       norm of the rows of A. It should be used as a left preconditioner. 
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

function PL = mpi_diagonal_preconditioner (A, Comm)
  Adim = size(A);
  m    = Adim(1);
  row  = zeros(1,m);
  for i=1:m
    row      = MPI_Allreduce(A(i,:),'mpi_sum',Comm);
    row_norm = sum(row.^2);
    D(i)     = sqrt(row_norm);
  end
  PL = sparse( diag( D ) );
end