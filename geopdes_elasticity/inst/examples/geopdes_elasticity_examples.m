% GEOPDES_ELASTICITY_EXAMPLES: Run some simple examples on how to use the geopdes_elasticity package.
%
% Copyright (C) 2006-2009, Thomas Treichl <treichl@users.sourceforge.net>
% Copyright (C) 2010, 2011 Rafael Vazquez
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
            ' 2D Plane strain examples. \n \n', ...
            '   (1) Unit square with Dirichlet boundary conditions. \n \n', ...
            '   (2) Unit square with mixed boundary conditions. \n \n', ...
            '   (3) Quarter of a ring. \n \n', ...
            '   (4) Plate with a hole. \n \n', ...
            '\n', ...
            ' 3D Linear elasticity examples. \n \n', ...
            '   (5) Horseshoe (thanks to T.J.R. Hughes and his group). \n \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

  if (iopt == 1)
    clc;
    fprintf (1, 'You can have a look at the source file: EX_PLANE_STRAIN_SQUARE \n \n');
    ex_plane_strain_square;
    input ('Press <Enter> to continue: ');
  elseif (iopt == 2)
    clc;
    fprintf (1, 'You can have a look at the source file: EX_PLANE_STRAIN_SQUARE_MIXED_BC \n \n');
    ex_plane_strain_square_mixed_bc;
    input ('Press <Enter> to continue: ');
  elseif (iopt == 3)
    clc;
    fprintf (1, 'You can have a look at the source file: EX_PLANE_STRAIN_RING \n \n');
    ex_plane_strain_ring;
    input ('Press <Enter> to continue: ');
  elseif (iopt == 4)
    clc;
    fprintf (1, 'You can have a look at the source file: EX_PLANE_STRAIN_PLATE \n \n');
    ex_plane_strain_plate;
    input ('Press <Enter> to continue: ');
  elseif (iopt == 5)
    clc;
    fprintf (1, 'You can have a look at the source file: EX_LIN_ELAST_HORSESHOE \n \n');
    ex_lin_elast_horseshoe;
    input ('Press <Enter> to continue: ');
  end
end %# while (iopt > 0)

end