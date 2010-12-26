% MPI_ALLREDUCE: emulation of the equally named mpi function
%
% result = MPI_Allreduce(value, op, comm)
%
% INPUT:
%
% value:      the value to send,
% op:         the merge operation. The supported operations are
%             'mpi_max', 'mpi_min', 'mpi_sum' and 'mpi_prod'. 
%             The operations are case insensitive.
% comm:       mpi communication group
%
% OUTPUT:
%
% result:     the result of the operation applied to the the values
%             supplied by all processors
%
% NOTES:
%
% This function does not handle errors nor checks for them. If a
% communication fails the error will appear when applying the operation.
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



function result=MPI_Allreduce ( myvalue, op, comm)

TAG    = 716253; % random number, fixed to leave the interface as
		 % similar as possible to that of mpi API

cpu_N  = MPI_Comm_size ( comm );

if cpu_N == 1
  result=myvalue;
  return;
end

cpu_ID = MPI_Comm_rank ( comm );

nel    = prod(size (myvalue));
shape  = size (myvalue);
recv   = zeros(cpu_N , nel);

for i=1:cpu_ID
  [temp info] = MPI_Recv(i-1,TAG,comm);
  recv(i,:) = temp(:);
end

dest = setdiff (0:cpu_N-1, cpu_ID);
MPI_Send(myvalue, dest, TAG, comm);
recv(cpu_ID+1,:)=myvalue(:);

for i=cpu_ID+2:cpu_N
  temp = MPI_Recv(i-1,TAG,comm);
  recv(i,:) = temp(:);
end

switch lower(op)
  case 'mpi_max'
    result = max (recv);
  case 'mpi_min'
    result = min (recv);
  case 'mpi_sum'
    result = sum (recv);
  case 'mpi_prod'
    result = prod (recv);
  otherwise
    error('uniplemented MPI_OPERATION: %s', op);
end
result=reshape(result, shape);
end

%!test
%!
%! MPI_Init();
%!
%! comm   = MPI_Comm_Load('ciao');
%! cpu_N  = MPI_Comm_size(comm);
%! cpu_ID = MPI_Comm_rank(comm);
%!
%! a = [1 2 3]*(cpu_ID+1);
%! b = MPI_Allreduce(a,'mpi_sum' ,comm);
%! sleep(0.3*cpu_ID);
%! printf('Rank: %d; Sum: ',cpu_ID)
%! disp(b)
%! assert (b== [1 2 3]*sum(1:cpu_N));
%!
%! a = [1 .5 .33];
%! b =  MPI_Allreduce(a,'mpi_prod' ,comm);
%! sleep(0.3*cpu_ID);
%! printf('Rank: %d; Prod: ',cpu_ID)
%! disp(b)
%! assert (b==a.^cpu_N);
%! 
%! a = [1 2 3]*(cpu_ID+1);
%! b =  MPI_Allreduce(a,'mpi_max' ,comm);
%! sleep(0.3*cpu_ID);
%! printf('Rank: %d; Max: ',cpu_ID)
%! disp(b)
%! assert (b==[1 2 3]*cpu_N);
%! 
%! a = [1 2 3]*(cpu_ID+1);
%! b =  MPI_Allreduce(a,'mpi_min' ,comm);
%! sleep(0.3*cpu_ID);
%! printf('Rank: %d; Min: ',cpu_ID)
%! disp(b)
%! assert (b==[1 2 3]);
%! MPI_Finalize();
