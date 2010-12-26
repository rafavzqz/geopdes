% MPI_STATUS_MESSAGE: prints a message and the time of call, usefull for
%                     profiling in a parallel environment.
%                     Node 0 is assumed to have io capabilities.
%
%   mpi_timed_message( msg, comm)
%
% INPUT:
%
%  msg:  the message to be printed
%  comm: the mpi communication group
%    
%
% Copyright (C) 2010 Andrea Bressan
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

function mpi_timed_message(msg, Comm )

Rank=MPI_Comm_rank(Comm);

if Rank==0
  time_of_call = clock;
  time_string  = datestr(time_of_call);
  message=sprintf('At %s: %s', time_string, msg);
  disp(message);
end
end