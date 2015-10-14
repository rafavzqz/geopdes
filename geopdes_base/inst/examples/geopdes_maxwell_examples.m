% GEOPDES_MAXWELL_EXAMPLES: Run some simple examples on how to use the geopdes_maxwell package.
%
% Copyright (C) 2006-2009, Thomas Treichl <treichl@users.sourceforge.net>
% Copyright (C) 2010-2011, Rafael Vazquez
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

function [] = geopdes_maxwell_examples ()

iopt = 1; 
while (iopt > 0)
  clc;
  fprintf (1, ...
           ['GeoPDEs_Maxwell examples menu:\n', ...
            '------------------------------\n', ...
            '\n', ...
            '   (1) Examples in 2D: source problem. \n \n',...
            '   (2) Examples in 2D: eigenvalue problem. \n \n',...
            '   (3) Examples in 3D: source problem. \n \n',...
            '   (4) Examples in 3D: eigenvalue problem. \n \n', ...
            '   (5) Multipatch geometries: eigenvalue problem (geopdes_multipatch required). \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

  if (iopt == 2)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs_Maxwell examples menu: 2D eigenvalue problem\n', ...
                '----------------------------------------------------\n', ...
                '\n', ...
                'Examples with the direct formulation \n \n', ...
                '   (1) Unit square domain \n', ...
                '   (2) One eighth of a ring \n', ...
                '   (3) L-shaped domain \n', ...
                '   (4) Curved L-shaped domain \n \n', ...
                'Examples with the first mixed formulation \n \n', ...
                '   (5) Unit square domain \n', ...
                '   (6) L-shaped domain \n \n', ...
                'Examples with the second mixed formulation \n \n', ...
                '   (7) Unit square domain \n', ...
                '   (8) L-shaped domain \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_square.m \n \n');
          ex_maxwell_eig_square;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_ring_1eighth.m \n \n');
          ex_maxwell_eig_ring_1eighth;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_Lshaped.m \n \n');
          ex_maxwell_eig_Lshaped;
          input ('Press <Enter> to continue: ');
         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_curvedL.m \n \n');
          ex_maxwell_eig_curvedL;
          input ('Press <Enter> to continue: ');
         case 5
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_square.m \n \n');
          ex_maxwell_eig_mixed1_square;
          input ('Press <Enter> to continue: ');
         case 6
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_Lshaped.m \n \n');
          ex_maxwell_eig_mixed1_Lshaped;
          input ('Press <Enter> to continue: ');
         case 7
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed2_square.m \n \n');
          ex_maxwell_eig_mixed2_square;
          input ('Press <Enter> to continue: ');
         case 8
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed2_Lshaped.m \n \n');
          ex_maxwell_eig_mixed2_Lshaped;
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)

  elseif (iopt == 4)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs_Maxwell examples menu: 3D eigenvalue problem\n', ...
                '----------------------------------------------------\n', ...
                '\n', ...
                'Examples with the direct formulation \n \n', ...
                '   (1) Unit cube domain \n', ...
                '   (2) Thick L-shaped domain \n \n', ...
                'Examples with the first mixed formulation \n \n', ...
                '   (3) Unit cube domain \n', ...
                '   (4) Thick L-shaped domain \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_cube.m \n \n');
          ex_maxwell_eig_cube;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_thick_L.m \n \n');
          ex_maxwell_eig_thick_L;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_cube.m \n \n');
          ex_maxwell_eig_mixed1_cube;
          input ('Press <Enter> to continue: ');
         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_thick_L.m \n \n');
          ex_maxwell_eig_mixed1_thick_L;
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)

  elseif (iopt == 1)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
         ['GeoPDEs_Maxwell examples menu: 2D source problem\n', ...
          '------------------------------------------------\n', ...
          '\n', ...
          ' 2D Maxwell source problem solved with the direct formulation. \n \n', ...
          '   (1) Unit square domain \n \n', ...
          '   (2) L-shaped domain \n \n', ...
          '   (3) One quarter of a ring \n \n', ...
          '   (4) Three quarters of a circle \n \n']);

     iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_src_square.m \n \n');
          ex_maxwell_src_square;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_src_Lshaped.m \n \n');
          ex_maxwell_src_Lshaped;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_src_ring.m \n \n');
          ex_maxwell_src_ring;
          input ('Press <Enter> to continue: ');
         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_src_pacman.m \n \n');
          ex_maxwell_src_pacman;
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)

  elseif (iopt == 3)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
          ['GeoPDEs_Maxwell examples menu: 3D source problem\n', ...
           '------------------------------------------------\n', ...
           '\n', ...
           '   3D Maxwell source problem solved with the direct formulation. \n \n', ...
          '   (1) Unit cube domain \n \n', ...
          '   (2) One quarter of a thick ring \n \n', ...
          '   (3) Three quarters of a cylinder \n \n']);

     iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_src_cube.m \n \n');
          ex_maxwell_src_cube;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_src_thick_ring.m \n \n');
          ex_maxwell_src_thick_ring;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_src_camembert.m \n \n');
          ex_maxwell_src_camembert;
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)

  elseif (iopt == 5)
    if (~exist('mp_geo_load'))
      fprintf(1, 'Unable to find the file ''mp_geo_load''.\n');
      fprintf(1, 'Be sure to install the latest version of geopdes_multipatch to run these examples\n\n');
      iopt2 = -1;
      input ('Press <Enter> to continue: ');
    else
      iopt2 = 1;
    end
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
        ['GeoPDEs_maxwell examples menu: eigenvalue problem in multipatch geometries\n', ...
         '-------------------------------------------------------------------------- \n', ...
         '\n', ...
         '2D examples \n \n', ...
         '   (1) L-shaped domain, defined with 3 patches. \n',...
         '   (2) L-shaped domain, defined with 3 patches. Mixed formulation. \n \n',...
         '3D examples \n \n', ...
         '   (3) Unit cube, defined with 2 patches. \n',...
         '   (4) Unit cube, defined with 2 patches. Mixed formulation.\n\n',...
         '   (5) Thick L-shaped domain, defined with 3 patches. \n',...
         '   (6) Thick L-shaped domain, defined with 3 patches. Mixed formulation.\n\n',...
         '   (7) Fichera''s corner, defined with 7 patches. \n',...
         '   (8) Fichera''s corner, defined with 7 patches. Mixed formulation.\n \n']);

      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_Lshaped_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_Lshaped_mp;
          input ('Press <Enter> to continue: ');

         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_Lshaped_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_mixed1_Lshaped_mp;
          input ('Press <Enter> to continue: ');

         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_cube_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_cube_mp;
          input ('Press <Enter> to continue: ');

         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_cube_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_mixed1_cube_mp;
          input ('Press <Enter> to continue: ');

         case 5
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_thick_L_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_thick_L_mp;
          input ('Press <Enter> to continue: ');

         case 6
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_thick_L_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_mixed1_thick_L_mp;
          input ('Press <Enter> to continue: ');

         case 7
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_fichera_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_fichera_mp;
          input ('Press <Enter> to continue: ');

         case 8
          clc;
          fprintf (1, 'You can have a look at the source file: ex_maxwell_eig_mixed1_fichera_mp.m \n \n');
          fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
          ex_maxwell_eig_mixed1_fichera_mp;
          input ('Press <Enter> to continue: ');
        end %# switch
      end %# if (~isempty)
    end %# while iopt2>0

  end

end %# while (iopt > 0)
end

