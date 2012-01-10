% SP_TO_VTK_STRESS: Export the components of the stress tensor to VTK format for plotting.
%
%  sp_to_vtk_stress (u, space, geometry, npts, lambda, mu, filename)
%  sp_to_vtk_stress (u, space, geometry, pts, lambda, mu, filename)
%
% INPUT:
%     
%     u:          vector of dof weights
%     space:      object representing the space of discrete functions (see sp_vector_2d)
%     geometry:   geometry structure (see geo_load)
%     npts:       number of points along each parametric direction where to evaluate
%     pts:        cell array with the coordinates along each parametric direction of the points where to evaluate
%     lambda, mu: Lame' parameters
%     filename:   name of the output file. 
%
% OUTPUT:
%
%    none    
% 
% Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011, 2012 Rafael Vazquez
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

function sp_to_vtk_stress (u, space, geometry, npts, lambda, mu, filename)

  [stress, pts] = sp_eval_stress (u, space, geometry, npts, lambda, mu);

  if (length (size (pts)) == 4)
    dim = 3;
  else
    dim = 2;
  end

  str1 = cat (2,'<?xml version="1.0"?> \n', ...
'<VTKFile type="StructuredGrid" version="0.1"> \n', ...
'<StructuredGrid WholeExtent="0 %d 0 %d 0 %d"> \n', ...
'<Piece Extent="0 %d 0 %d 0 %d"> \n', ...
'<PointData Scalars="Stress components">\n');

  str2 = ...
'<DataArray type="Float32" Name="%s" format="ascii" NumberOfComponents="1"> \n';

  str3 = '</DataArray> \n';

  str4 = cat (2,'</PointData> \n', ...
'<Points> \n', ...
'<DataArray type="Float32" NumberOfComponents="3"> \n');

  str5 = cat (2, '\n', ...
'</DataArray>\n', ...
'</Points> \n', ...
'</Piece> \n', ...
'</StructuredGrid> \n', ...
'</VTKFile> \n');

  if (dim == 3)
    npts = size (squeeze (pts(1,:,:,:)));
  else
    npts = size (squeeze (pts(1,:,:)));
    npts(3) = 1;
    pts(3,:,:) = 0;
  end

  if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.vts'))
    filename = cat (2, filename, '.vts');
  end

  fid = fopen (filename, 'w');
  if (fid < 0)
    error ('sp_to_vtk_stress: could not open file %s', filename);
  end

  fprintf (fid, str1, ...
           npts(1)-1, npts(2)-1, npts(3)-1, ...
           npts(1)-1, npts(2)-1, npts(3)-1);

  if (space.ncomp == 2)
    fprintf (fid, str2, 'sigma_xx');
    fprintf (fid, '%g ', reshape (squeeze (stress(1,1,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str2, 'sigma_xy');
    fprintf (fid, '%g ', reshape (squeeze (stress(1,2,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str2, 'sigma_yy');
    fprintf (fid, '%g ', reshape (squeeze (stress(2,2,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str4);
    fprintf (fid, '%g ', pts(:));
    fprintf (fid, str5);

    fclose (fid);
  elseif (space.ncomp == 3)
    fprintf (fid, str2, 'sigma_xx');
    fprintf (fid, '%g ', reshape (squeeze (stress(1,1,:,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str2, 'sigma_xy');
    fprintf (fid, '%g ', reshape (squeeze (stress(1,2,:,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str2, 'sigma_xz');
    fprintf (fid, '%g ', reshape (squeeze (stress(1,3,:,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str2, 'sigma_yy');
    fprintf (fid, '%g ', reshape (squeeze (stress(2,2,:,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str2, 'sigma_yz');
    fprintf (fid, '%g ', reshape (squeeze (stress(2,3,:,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str2, 'sigma_zz');
    fprintf (fid, '%g ', reshape (squeeze (stress(3,3,:,:,:)), [], 1));
    fprintf (fid, str3);
    fprintf (fid, str4);
    fprintf (fid, '%g ', pts(:));
    fprintf (fid, str5);

    fclose (fid);

  end

end
