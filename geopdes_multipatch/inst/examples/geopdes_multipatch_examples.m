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
            '   (1) Poisson problem. \n \n',...
            '   (2) Linear elasticity problems (geopdes_elasticity must be installed) \n \n',...
            '   (3) Stokes problem (geopdes_fluid must be installed) \n \n',...
            '   (4) Maxwell eigenvalue problem (geopdes_maxwell must be installed). \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

  if (iopt == 1)
    iopt2 = 1;
    while (iopt2 > 0)
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

      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
      switch iopt2
       case 1
        clc;
        fprintf (1, 'You can have a look at the source file: EX_LAPLACE_BSP_LSHAPED_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_laplace_bsp_Lshaped_mp;
        input ('Press <Enter> to continue: ');

       case 2
        clc;
        fprintf (1, 'You can have a look at the source file: EX_LAPLACE_BSP_CUBE_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_laplace_bsp_cube_mp;
        input ('Press <Enter> to continue: ');

       case 3
        clc;
        fprintf (1, 'You can have a look at the source file: EX_LAPLACE_BSP_THICK_L_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_laplace_bsp_thick_L_mp;
        input ('Press <Enter> to continue: ');
      end %switch
    end %while iopt2>0

  elseif (iopt == 2)
    if (~exist('ex_plane_strain_Lshaped_mp'))
      fprintf(1, 'You must install geopdes_elasticity to run the examples\n\n');
      iopt2 = -1;
      input ('Press <Enter> to continue: ');
    else
      iopt2 = 1;
    end
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
        ['GeoPDEs_multipatch examples menu: Linear elasticity problems \n', ...
         '------------------------------------------------------------\n', ...
         '\n', ...
         '2D examples: plane strain \n \n', ...
         '   (1) L-shaped domain, defined with 3 patches. \n \n',...
         '3D examples \n \n', ...
         '   (2) Unit cube, defined with 2 patches. \n \n']);

      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      switch iopt2
       case 1
        clc;
        fprintf (1, 'You can have a look at the source file: EX_PLANE_STRAIN_LSHAPED_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_plane_strain_Lshaped_mp;
        input ('Press <Enter> to continue: ');

       case 2
        clc;
        fprintf (1, 'You can have a look at the source file: EX_LIN_ELAST_CUBE_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_lin_elast_cube_mp;
        input ('Press <Enter> to continue: ');
      end
    end %while iopt2>0

  elseif (iopt == 3)
    if (~exist('ex_stokes_mp'))
      fprintf(1, 'GUARDA IL FILE DEGLI ESEMPI \n\n');
      fprintf(1, 'You must install geopdes_fluid to run the examples\n\n');
      iopt2 = -1;
      input ('Press <Enter> to continue: ');
    else
      iopt2 = 1;
    end
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
        ['GeoPDEs_multipatch examples menu: Stokes problem \n', ...
         '------------------------------------------------ \n', ...
         '\n', ...
         '2D examples \n \n', ...
         '   (1) NON HO ANCORA FINITO \n \n',...
         '3D examples \n \n', ...
         '   (2) NON HO ANCORA FINITO \n \n']);

      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
    end %while iopt2>0


  elseif (iopt == 4)
    if (~exist('ex_maxwell_eig_cube_mp'))
      fprintf(1, 'You must install geopdes_maxwell to run the examples\n\n');
      iopt2 = -1;
      input ('Press <Enter> to continue: ');
    else
      iopt2 = 1;
    end
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
        ['GeoPDEs_multipatch examples menu: Maxwell eigenvalue problem \n', ...
         '------------------------------------------------------------ \n', ...
         '\n', ...
         '2D examples \n \n', ...
         '   (1) L-shaped domain, defined with 3 patches. \n',...
         '   (2) L-shaped domain, defined with 3 patches. Mixed formulation. \n \n',...
         '3D examples \n \n', ...
         '   (3) Unit cube, defined with 2 patches. \n',...
         '   (4) Unit cube, defined with 2 patches. Mixed formulation.\n',...
         '   (5) Thick L-shaped domain, defined with 3 patches. \n',...
         '   (6) Thick L-shaped domain, defined with 3 patches. Mixed formulation.\n',...
         '   (7) Fichera corner, defined with 7 patches. \n',...
         '   (8) Fichera corner, defined with 7 patches. Mixed formulation.\n \n']);

      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');


      switch iopt2
       case 1
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_LSHAPED_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_Lshaped_mp;
        input ('Press <Enter> to continue: ');

       case 2
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_LSHAPED_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_mixed1_Lshaped_mp;
        input ('Press <Enter> to continue: ');

       case 3
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_CUBE_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_cube_mp;
        input ('Press <Enter> to continue: ');

       case 4
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_CUBE_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_mixed1_cube_mp;
        input ('Press <Enter> to continue: ');

       case 5
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_THICK_L_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_thick_L_mp;
        input ('Press <Enter> to continue: ');

       case 6
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_THICK_L_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_mixed1_thick_L_mp;
        input ('Press <Enter> to continue: ');

       case 7
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_FICHERA_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_fichera_mp;
        input ('Press <Enter> to continue: ');

       case 8
        clc;
        fprintf (1, 'You can have a look at the source file: EX_MAXWELL_EIG_MIXED1_FICHERA_MP \n \n');
        fprintf (1, 'You may also modify the file to solve in the same geometry with rotated patches\n \n');
        ex_maxwell_eig_mixed1_fichera_mp;
        input ('Press <Enter> to continue: ');

    end %while iopt2>0

  end

end %# while (iopt > 0)

end

