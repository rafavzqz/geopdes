% MPI_WORK_INIT: initialize the mpi environment
%
%  [Comm_group]=mpi_work_init()
%
% INPUT:
%
% OUTPUT:
%
%	comm_group:	mpi communication group necessary to call any mpi function
%
% Comm_group: a structure containing the MPI_Comm_Group data
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

function  Comm_group=mpi_work_init()
	MPI_SUCCESS=0;
	if ~MPI_Initialized()
		error_code=MPI_Init();
		if error_code ~=  MPI_SUCCESS
			error('MPI could not be initialized, check your configuration');
		end
	end
	Comm_group	=MPI_Comm_Load('All nodes');
end
