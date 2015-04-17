% MSH_TO_VTK: Export to VTK format for plotting.
%
%  MSH_to_vtk (pts, values, filename, fieldname)
%
% INPUT:
%
%     pts:       points at which the field was computed
%     values:    values of the field at the selected point
%     filename:  name of the output file
%     fieldname: how to name the saved variable in the vtk file
%
% OUTPUT:
%
%    a vtk structured mesh file named <filename> is produced 
% 
% Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011 Rafael Vazquez
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

function msh_to_vtk (pts, values, filename, fieldname)

  ndim = numel (size (pts)) - 1;
  rdim = size (pts, 1);
  
  str1 = cat (2,'<?xml version="1.0"?> \n', ...
'<VTKFile type="StructuredGrid" version="0.1"> \n', ...
'<StructuredGrid WholeExtent="0 %d 0 %d 0 %d"> \n', ...
'<Piece Extent="0 %d 0 %d 0 %d"> \n', ...
'<PointData %s="%s">\n', ...
'<DataArray type="Float32" Name="%s" format="ascii" NumberOfComponents="%d"> \n');

  str2 = cat (2,'</DataArray> \n', ...
'</PointData> \n', ...
'<Points> \n', ...
'<DataArray type="Float32" NumberOfComponents="3"> \n');

  str3 = cat (2, '\n', ...
'</DataArray>\n', ...
'</Points> \n', ...
'</Piece> \n', ...
'</StructuredGrid> \n', ...
'</VTKFile> \n');

% Even for 2D data, everything is saved in 3D 
% ndims (or size) do not work properly for the 1D case. I remove singleton
% dimensions using this trick
  if (sum (size (values) > 1) == ndim)
    fieldclass = 'Scalars';
    ncomp = 1;
  else
    fieldclass = 'Vectors';
    ncomp = 3;
  end

  size_pts = size (pts);
  npts = size_pts (2:end);
  
  if (ndim < 3)
    npts (ndim+1:3) = 1;
  end
  if (rdim < 3)
    pts(ndim+1:3,:,:) = 0;
  end
  if (ncomp == 3 && rdim < 3)
    values(rdim+1:3,:,:) = 0;
  end

  if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.vts'))
    filename = cat (2, filename, '.vts');
  end

  fid = fopen (filename, 'w');
  if (fid < 0)
    error ('msh_to_vtk: could not open file %s', filename);
  end

  fprintf (fid, str1, ...
           npts(1)-1, npts(2)-1, npts(3)-1, ...
           npts(1)-1, npts(2)-1, npts(3)-1,...
           fieldclass, fieldname, fieldname, ncomp);
  
  fprintf (fid, '%g ', values(:));
  fprintf (fid, str2);
  fprintf (fid, '%g ', pts(:));
  fprintf (fid, str3);

  fclose (fid);

end
