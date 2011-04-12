% SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = sp_drchlt_l2_proj (sp, msh, h, sides)
%
% INPUT:
%
%  sp:    structure representing the space of trial functions (see sp_bspline_2d_phys)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%  h:     function handle to compute the Dirichlet condition
%  sides: boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

function [u, dofs] = mpi_sp_drchlt_l2_proj (mpi_comm, mpi_prec, mpi_solv, sp, msh, h, sides)

  dofs = unique ([sp.boundary(sides).dofs]);
  dofs = MPI_Allmerge( dofs, mpi_comm);
  
  M    = spalloc (sp.ndof, sp.ndof, sp.ndof);
  rhs  = spalloc (sp.ndof, 1, numel (dofs));
  
  for iside = intersect (sides, msh.boundary_list)
    
    msh_bnd = msh.boundary(iside);
    sp_bnd  = sp.boundary(iside);

    if (size (msh_bnd.geo_map, 1) == 2)
      [x, y] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                     squeeze (msh_bnd.geo_map(2,:,:)));

      hval = reshape (h (x, y, iside), sp.ncomp, msh_bnd.nqn, msh_bnd.nel);

    elseif (size (msh_bnd.geo_map, 1) == 3)
      [x, y, z] = deal (squeeze (msh_bnd.geo_map(1,:,:)), ...
                        squeeze (msh_bnd.geo_map(2,:,:)), ...
                        squeeze (msh_bnd.geo_map(3,:,:)));

      hval = reshape (h (x, y, z, iside), sp_bnd.ncomp, msh_bnd.nqn, msh_bnd.nel);
    end

    M_side = op_u_v (sp_bnd, sp_bnd, msh_bnd, ones (size (x)));
    M(sp_bnd.dofs, sp_bnd.dofs) = M(sp_bnd.dofs, sp_bnd.dofs) + M_side;

    rhs_side = op_f_v (sp_bnd, msh_bnd, hval);
    rhs(sp_bnd.dofs) = rhs(sp_bnd.dofs) + rhs_side;

  end
  rhs = MPI_Allreduce(rhs , 'mpi_sum', mpi_comm);
  
  preconditioner  =  mpi_prec( M(dofs, dofs), mpi_comm );
  u               =  mpi_solv( M(dofs,dofs), rhs(dofs), preconditioner, mpi_comm);
end
