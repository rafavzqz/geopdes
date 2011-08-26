% GEOPDES_FLUID_EXAMPLES: Run some simple examples on how to use the geopdes_fluid package.
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

function [] = geopdes_fluid_examples ()

iopt = 1; 
while (iopt > 0)
  clc;
  fprintf (1, ...
           ['GeoPDEs_Fluid examples menu:\n', ...
            '----------------------------\n', ...
            '\n', ...
            '   (1) Stokes flow: 2D examples. \n \n', ...
            '   (2) Stokes flow: 3D examples. \n \n', ...
            '   (3) Stokes flow in multipatch geometries (require geopdes_multipatch). \n \n']);
%            '3D Stokes flow examples. \n \n', ...
%            '   (9) 3D driven cavity discretized with generalized Taylor-Hood elements. \n',...
%            '   (10) 3D driven cavity discretized with the subgrid method. \n',...
%            '   (11) Twisted pipe discretized with generalized Taylor-Hood elements. \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

  if (iopt == 1)
    iopt2 = 1;
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs_Fluid examples menu: 2D Stokes problem\n', ...
                '----------------------------------------------\n', ...
                '\n', ...
                '   (1) Unit square, discretized with generalized Taylor-Hood elements.\n', ...
                '   (2) Unit square, discretized with generalized Raviart-Thomas elements.\n', ...
                '   (3) Unit square, discretized with generalized Nedelec elements.\n', ...
                '   (4) Unit square, discretized with the subgrid method.\n \n', ...
                '   (5) Quarter of a ring, discretized with generalized Taylor-Hood elements.\n', ...
                '   (6) Quarter of a ring, discretized with generalized Raviart-Thomas elements.\n\n', ...
                '   (7) Driven cavity discretized with generalized Taylor-Hood elements. \n',...
                '   (8) Driven cavity discretized with generalized Raviart-Thomas elements. \n\n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
      
      if (~isempty (iopt2) && iopt2 > 0 && iopt2 < 9)
        switch iopt2
         case 1
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_square_th.m \n \n');
          ex_stokes_square_th;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_square_rt.m \n \n');
          ex_stokes_square_rt;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_square_ndl.m \n \n');
          ex_stokes_square_ndl;
          input ('Press <Enter> to continue: ');
         case 4
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_square_sg.m \n \n');
          ex_stokes_square_sg;
          input ('Press <Enter> to continue: ');
         case 5
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_annulus_th.m \n \n');
          ex_stokes_annulus_th;
          input ('Press <Enter> to continue: ');
         case 6
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_annulus_rt.m \n \n');
          ex_stokes_annulus_rt;
          input ('Press <Enter> to continue: ');
         case 7
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_article_ijnmf_th.m \n \n');
          ex_stokes_article_ijnmf_th;
          input ('Press <Enter> to continue: ');
         case 8
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_article_ijnmf_rt.m \n \n');
          ex_stokes_article_ijnmf_rt;
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)

  elseif (iopt == 2)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs_Fluid examples menu: 3D problems\n', ...
                '----------------------------------------\n', ...
                '\n', ...
                '   (1) Symmetric driven cavity problem. Generalized Taylor-Hood elements.\n \n',...
                '   (2) Symmetric driven cavity problem. Subgrid method.\n \n',...
                '   (3) Twisted pipe. Generalized Taylor-Hood elements. \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
      
      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_driven_cavity_3d_th.m \n \n');
          ex_stokes_driven_cavity_3d_th;
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_driven_cavity_3d_sg.m \n \n');
          ex_stokes_driven_cavity_3d_sg;
          input ('Press <Enter> to continue: ');
         case 3
          clc;
          fprintf (1, 'You can have a look at the source file: ex_stokes_twisted_pipe.m \n \n');
          ex_stokes_twisted_pipe;
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
        
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)


  elseif (iopt == 3)
    clc;
    if (~exist('mp_interface_vector_2d'))
      fprintf(1, 'Unable to find the file ''mp_interface_vector_2d''.\n');
      fprintf(1, 'Be sure to install the latest version of geopdes_multipatch to run these examples\n\n');
      input ('Press <Enter> to continue: ');
    else
      iopt2 = 1; 
      while (iopt2 > 0)
        clc;
        fprintf (1, ...
               ['GeoPDEs_Fluid examples menu: multipatch problems\n', ...
                '------------------------------------------------\n', ...
                '\n', ...
                '   (1) 2D bifurcation problem. Generalized Taylor-Hood elements.\n \n',...
                '   (2) 3D multipatch driven cavity. Generalized Taylor-Hood elements. \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
      
        if (~isempty (iopt2))
          switch iopt2
           case 1
            clc;
            fprintf (1, 'You can have a look at the source file: ex_stokes_bifurcation_2d_mp.m \n \n');
            ex_stokes_bifurcation_2d_mp;
            input ('Press <Enter> to continue: ');
           case 2
            clc;
            fprintf (1, 'You can have a look at the source file: ex_stokes_driven_cavity_3d_mp.m \n \n');
            ex_stokes_driven_cavity_3d_mp;
            input ('Press <Enter> to continue: ');
          end %# switch (iopt2)
        end %# if (~isempty)
      end %# while (iopt2 > 0)
    end %# if (~exist())
  end
  
end %# while (iopt > 0)
end
