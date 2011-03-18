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
            '   (4) Examples in 3D: eigenvalue problem. \n \n']);

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
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_SQUARE \n \n');
          ex_maxwell_eig_square;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_RING_1EIGHTH \n \n');
          ex_maxwell_eig_ring_1eighth;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_LSHAPED \n \n');
          ex_maxwell_eig_Lshaped;
          input ('Press <Enter> to continue: ');
         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_CURVEDL \n \n');
          ex_maxwell_eig_curvedL;
          input ('Press <Enter> to continue: ');
         case 5
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_SQUARE \n \n');
          ex_maxwell_eig_mixed1_square;
          input ('Press <Enter> to continue: ');
         case 6
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_LSHAPED \n \n');
          ex_maxwell_eig_mixed1_Lshaped;
          input ('Press <Enter> to continue: ');
         case 7
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED2_SQUARE \n \n');
          ex_maxwell_eig_mixed2_square;
          input ('Press <Enter> to continue: ');
         case 8
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED2_LSHAPED \n \n');
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
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_CUBE \n \n');
          ex_maxwell_eig_cube;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_THICK_L \n \n');
          ex_maxwell_eig_thick_L;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_CUBE \n \n');
          ex_maxwell_eig_mixed1_cube;
          input ('Press <Enter> to continue: ');
         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_THICK_L \n \n');
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
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_SRC_SQUARE \n \n');
          ex_maxwell_src_square;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_SRC_LSHAPED \n \n');
          ex_maxwell_src_Lshaped;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_SRC_RING \n \n');
          ex_maxwell_src_ring;
          input ('Press <Enter> to continue: ');
         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_SRC_PACMAN \n \n');
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
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_SRC_CUBE \n \n');
          ex_maxwell_src_cube;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_SRC_THICK_RING \n \n');
          ex_maxwell_src_thick_ring;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: EX_MAXWELL_SRC_CAMEMBERT \n \n');
          ex_maxwell_src_camembert;
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)

  end

end %# while (iopt > 0)
end

