% MPI_ALLMERGE: 
%
% result = MPI_Allmerge(set, comm)
%
% INPUT:
%
% value:      the value to send,
% comm:       mpi communication group
%
% OUTPUT:
%
% result:     the union of the sets supplied by all processors
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



function result=MPI_Allmerge ( myvalue, comm)

TAG    = 762346; % can be anything

cpu_N  = MPI_Comm_size ( comm );
cpu_ID = MPI_Comm_rank ( comm );

result=myvalue;

for i=1:cpu_ID
  temp   = MPI_Recv(i-1,TAG,comm);
  result = unique([ temp(:); result(:)]);
end

dest = setdiff (0:cpu_N-1, cpu_ID);
MPI_Send(myvalue, dest, TAG, comm);

for i=cpu_ID+2:cpu_N
  temp = MPI_Recv(i-1,TAG,comm);
  result = unique([ temp(:); result(:)]);
end
result = reshape(result,[1,length(result)]);
end

%!test
%! if ~MPI_Initialized()
%!    MPI_Init();
%! end
%! comm   = MPI_Comm_Load('ciao');
%! cpu_N  = MPI_Comm_size(comm);
%! cpu_ID = MPI_Comm_rank(comm);
%! t = 2*cpu_ID;
%! a = [t:t+1];
%! b = MPI_Allmerge(a,comm);
%! sleep(0.3*cpu_ID);
%! disp(b)
%! assert (b==0:2*cpu_N-1);
%! MPI_Finalize();
