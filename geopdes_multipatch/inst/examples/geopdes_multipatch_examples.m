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
           ['GeoPDEs_multipatch examples menu:\n', ...
            '---------------------------------\n', ...
            '\n', ...
            '   (1) Examples in 2D: Poisson problem. \n \n',...
            '   (2) Examples in 3D: Poisson problem. \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

  if (iopt == 1)
    clc;
    fprintf (1, ...
            ['GeoPDEs_multipatch examples menu: 2D examples\n', ...
             '---------------------------------------------\n', ...
             '\n', ...
             '   2D Poisson problem in a multipatch domain. \n \n', ...
             '   You can take a look at the script file: examples/ex_bspline_laplace_2d_mp.m \n \n', ...
             '   You can solve with one of the following data files: \n \n', ...
             '   - test_Lshaped_mp \n \n']);
    fprintf ('Please write the name of the data file, including quotation marks: \n');
    data_file = input ('[test_Lshaped_mp] ', 's');
    if (~isempty (data_file))
      eval(data_file);
    else 
      eval('test_Lshaped_mp');
    end
    vexa = do_example (1);
    eval (vexa);
fprintf (1, ...
        ['\nYou may modify the data file to solve in the same geometry with rotated patches\n']);
    input ('Press <Enter> to continue: ');

  elseif (iopt == 2)
    clc;
    fprintf (1, ...
            ['GeoPDEs_multipatch examples menu: 3D examples\n', ...
             '---------------------------------------------\n', ...
             '\n', ...
             '   3D Poisson problem in a multipatch domain. \n \n', ...
             '   You can take a look at the script file: examples/ex_bspline_laplace_3d_mp.m \n \n', ...
             '   You can solve with one of the following data files: \n \n', ...
             '   - test_cube_mp \n',...
             '   - test_thick_Lshaped_mp \n \n']);
    fprintf ('Please write the name of the data file, including quotation marks: \n');
    data_file = input ('[test_thick_Lshaped_mp] ', 's');
    if (~isempty (data_file))
      eval(data_file);
    else 
      eval('test_thick_Lshaped_mp');
    end
    vexa = do_example (2);
    eval (vexa);
    input ('Press <Enter> to continue: ');

  end

end %# while (iopt > 0)

end

function vexa = do_example (number)

switch (number)
 case 1
  filename = 'ex_bspline_laplace_2d_mp.m';
 case 2
  filename = 'ex_bspline_laplace_3d_mp.m';
end

fid = fopen (filename);
vexa = '';
l = fgets (fid);
while (l ~= -1)
  vexa = cat (2, vexa, l);
  l = fgets (fid);
end
end

