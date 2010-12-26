% GEOPDES_ELASTICITY_EXAMPLES: Run some simple examples on how to use the geopdes_elasticity package.
%
% Copyright (C) 2006-2009, Thomas Treichl <treichl@users.sourceforge.net>
% Copyright (C) 2010 Rafael Vazquez
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

function [] = geopdes_elasticity_examples ()

iopt = 1; 
while (iopt > 0)
  clc;
  fprintf (1, ...
           ['GeoPDEs_Elasticity examples menu:\n', ...
            '---------------------------------\n', ...
            '\n', ...
            '   (1) Examples in 2D (plane strain). \n \n', ...
            '   (2) Examples in 3D. \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;


  if (iopt == 1)
    clc;
    fprintf (1, ...
            ['GeoPDEs_Elasticity examples menu: plane strain problems \n', ...
             '--------------------------------------------------------\n', ...
             '\n', ...
             '   You can solve with one of the following data files: \n \n', ...
             '   - test_plane_strain_square \n',...
             '   - test_plane_strain_square_mixed_bc \n',...
             '   - test_plane_strain_ring \n',...
             '   - test_plane_strain_plate \n \n']);
    fprintf ('Please write the name of the data file: \n');
    data_file = input ('[test_plane_strain_plate] ', 's');
    if (isempty (data_file))
      data_file = 'test_plane_strain_plate';;
    end
    eval (data_file);
    fprintf (['\nRunning script ' '     ex_nurbs_plane_strain_2d.m', ...
              '\nwith the data file  %s \n', ...
              '\nYou can have a look at the source using the command:\n' ...
                       '   type ex_nurbs_plane_strain_2d \n\n'], data_file)
    ex_nurbs_plane_strain_2d
    input ('Press <Enter> to continue: ');

  elseif (iopt == 2)
    clc;
    fprintf (1, ...
            ['GeoPDEs_Elasticity examples menu: 3D linear elasticity problems\n', ...
             '----------------------------------------------------------------\n', ...
             '\n', ...
             '   You can solve with one of the following data files: \n \n', ...
             '   - test_linear_elasticity_horseshoe \n \n']);
    fprintf ('Please write the name of the data file: \n');
    data_file = input ('[test_linear_elasticity_horseshoe] ', 's');
    if (isempty (data_file))
      data_file = 'test_linear_elasticity_horseshoe';;
    end
    eval (data_file);
    fprintf (['\nRunning script ' '     ex_nurbs_linear_elasticity_3d.m', ...
              '\nwith the data file  %s \n', ...
              '\nYou can have a look at the source using the command:\n' ...
                       '   type ex_nurbs_linear_elasticity_3d \n\n'], data_file)
    ex_nurbs_linear_elasticity_3d
    input ('Press <Enter> to continue: ');
  end %# if (iopt)
end %# while (iopt > 0)
end

