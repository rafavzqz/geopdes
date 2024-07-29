% SP_COMPUTE_CPATCH_VECTOR: Compute the matrix for B-spline representation, and
%  the indices of C^1 functions that do not vanish on a given patch, for vector-valued functions.
%
% [Cpatch, Cpatch_cols] = sp_compute_Cpatch_vector (space, patch, ncomp)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_multipatch_C1)
%    patch:   index of the patch
%    ncomp:   number of components of the vector-valued space
%
% OUTPUT:
%    Cpatch:      coefficients of the linear combination of basis functions as standard B-splines,
%                   The matrix is block diagonal, with each block computed as sp_compute_Cpatch.
%    Cpatch_cols: indices of the C^1 basis functions that do not vanish on the patch
%                   (see also sp_get_functions_on_patch)
%
% Copyright (C) 2023 Rafael Vazquez
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

function [Cpatch, Cpatch_cols] = sp_compute_Cpatch_vector (space, iptc, rdim)

    [Cpatch, cols] = sp_compute_Cpatch (space, iptc);
    
    Cpatch = repmat ({Cpatch}, 1, rdim);
    Cpatch = blkdiag (Cpatch{:});
    
    Cpatch_cols = [];
    for icomp = 1:rdim
      Cpatch_cols = union (Cpatch_cols, (icomp-1)*space.ndof + cols);
    end

end