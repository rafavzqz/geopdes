% MP_SP_TO_VTK: Export multipatch results to VTK format for plotting.
%
%  mp_sp_to_vtk (u, space, geometry, gnum, npts, filename, fieldname)
%
% INPUT:
%     
%     u:         vector of dof weights
%     space:     structure representing the space of discrete functions (see sp_bspline_2d)
%     geometry:  geometry structure (see geo_load)
%     gnum:      array that relates the local numbering on each patch with the global one
%     npts:      number of points along each parametric direction where to evaluate
%     filename:  name of the output file. 
%     fieldname: how to name the saved variable in the vtk file
%
% OUTPUT:
%
%    none    
% 
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

function mp_sp_to_vtk (u, space, geometry, gnum, npts, filename, fieldname)

  str1 = cat (2,'<?xml version="1.0"?> \n', ...
'<VTKFile type="Collection" version="0.1"> \n', ...
'<Collection> \n');

  str2 = cat (2, '<DataSet part="%d" file="%s.vts"/> \n');

  str3 = cat (2, ...
'</Collection>\n', ...
'</VTKFile> \n');

  if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
    pvd_filename = cat (2, filename, '.pvd');
  else
    pvd_filename = filename;
    filename = filename (1:end-4);
  end

  fid = fopen (pvd_filename, 'w');
  if (fid < 0)
    error ('mp_sp_to_vtk: could not open file %s', pvd_filename);
  end

  fprintf (fid, str1);
  for iptc = 1:numel(geometry)
    filename_patch = cat (2, filename, '_', num2str (iptc));
    fprintf (fid, str2, iptc, filename_patch);
    sp_to_vtk (u(gnum{iptc}), space{iptc}, geometry(iptc), npts, ...
                           filename_patch, fieldname)
  end
  fprintf (fid, str3);

  fclose (fid);

end
