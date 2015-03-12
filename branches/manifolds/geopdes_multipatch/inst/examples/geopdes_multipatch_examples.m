% MULTIPATCH_EXAMPLES: Run some simple examples on how to use the geopdes_multipatch package.
%
% Copyright (C) 2006-2009, Thomas Treichl <treichl@users.sourceforge.net>
% Copyright (C) 2010-2011 Rafael Vazquez
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

function [] = geopdes_multipatch_examples ()

iopt = 1; 
while (iopt > 0)
  clc;
  fprintf (1, ...
        ['GeoPDEs_multipatch examples menu: Poisson problem \n', ...
         '------------------------------------------------- \n', ...
         '\n', ...
         '2D examples \n \n', ...
         '   (1) L-shaped domain, defined with 3 patches. \n \n',...
         '3D examples \n \n', ...
         '   (2) Unit cube, defined with 2 patches. \n \n',...
         '   (3) Thick L-shaped domain, defined with 3 patches. \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

    if (~isempty(iopt))
      switch iopt
       case 1
        clc;
        fprintf (1, 'You can have a look at the source file: ex_laplace_lshaped_mp \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_laplace_Lshaped_mp;
        input ('Press <Enter> to continue: ');

       case 2
        clc;
        fprintf (1, 'You can have a look at the source file: ex_laplace_cube_mp \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_laplace_cube_mp;
        input ('Press <Enter> to continue: ');

       case 3
        clc;
        fprintf (1, 'You can have a look at the source file: ex_laplace_thick_L_mp \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_laplace_thick_L_mp;
        input ('Press <Enter> to continue: ');
      end %switch
    end %if

end %# while (iopt > 0)

end

