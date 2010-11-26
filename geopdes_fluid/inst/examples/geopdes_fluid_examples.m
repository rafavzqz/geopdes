% GEOPDES_FLUID_EXAMPLES: Run some simple examples on how to use the geopdes_fluid package.
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

function [] = geopdes_fluid_examples ()

iopt = 1; 
while (iopt > 0)
  clc;
  fprintf (1, ...
           ['GeoPDEs_Fluid examples menu:\n', ...
            '----------------------------\n', ...
            '\n', ...
            '   (1) Examples in 2D. \n \n',...
            '   (2) Examples in 3D. \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;

  if (iopt == 1)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs_Fluid examples menu: 2D problems\n', ...
                '----------------------------------------\n', ...
                '\n', ...
                '   (1) Stokes problem on a square with homogeneous Dirichlet BCs.\n \n', ...
                '   (2) Stokes problem on a square with non-homogeneous Dirichlet BCs.\n \n', ...
                '   (3) Stokes problem on a NURBS geometry with homogeneous Dirichlet BCs.\n \n', ...
                '   (4) Symmetric driven cavity problem. \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
      
      if (~isempty (iopt2) && iopt2 > 0 && iopt2 < 5)
        switch iopt2
         case 1
          datafile = 'test_stokes_square';
          eval (datafile);
         case 2
          datafile = 'test_stokes_square_bc';
          eval (datafile);
         case 3
          datafile = 'test_stokes_annulus';
          eval (datafile);
         case 4
          datafile = 'test_stokes_symdrivcav';
          eval (datafile);
        end %# switch (iopt2)

        
        iopt3 = 1; 
        clc;
        fprintf (1, ...
                 ['GeoPDEs_Fluid examples menu: choose function spaces.\n', ...
                  '----------------------------------------------------\n', ...
                  '\n', ...
                  '   (1) TH function spaces.\n \n', ...
                  '   (2) NDL function spaces.\n \n', ...
                  '   (3) RT function spaces. \n \n']);
        
        iopt3 = input ('Please choose a number from above or press <Enter> to return: ');
        
        if (~isempty (iopt3))
          switch iopt3
           case 1
            if (~strcmpi (sp_type, 'th'))
              fun_space    = 'sp_bspline_th_2d_phys';
              der2         = false;
            end
            fprintf (['\nRunning script ' '     ex_bspline_stokes_2d.m', ...
                      '\nwith the data file  %s \n', ...
                      '\nYou can have a look at the source using the command:\n' ...
                      '   type ex_bspline_stokes_2d\n\n'], datafile)
            ex_bspline_stokes_2d;
           case 2
            if (~strcmpi (sp_type, 'ndl'))
              fun_space    = 'sp_bspline_ndl_2d_phys';
              der2         = true;
            end
            fprintf (['\nRunning script ' '     ex_bspline_stokes_2d.m', ...
                      '\nwith the data file  %s \n', ...
                      '\nYou can have a look at the source using the command:\n' ...
                      '   type ex_bspline_stokes_2d\n\n'], datafile)
            ex_bspline_stokes_2d;
           case 3
            if (~strcmpi (sp_type, 'rt'))
              fun_space    = 'sp_bspline_rt_2d_phys';
              der2         = true;
              if (~isempty (drchlt_sides) && any (degree ~= 3))
                degree     = [3 3];
                regularity = min (regularity, 2); 
                fprintf ('N.B. for RT elements with Dirichlet conditions only degree 3 may be used\n')
              end
            end
            fprintf (['\nRunning script ' '     ex_bspline_stokes_2d_rt.m', ...
                      '\nwith the data file  %s \n', ...
                      '\nYou can have a look at the source using the command:\n' ...
                      '   type ex_bspline_stokes_2d_rt\n\n'], datafile)
            ex_bspline_stokes_2d_rt;
          end %# switch (iopt3)
          input ('Press <Enter> to continue: ');
        end %# if (~isempty (iopt3))
      end %# while (iopt2 > 0)
    end %# if (~isempty (iopt2))

  elseif (iopt == 2)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs_Fluid examples menu: 2D problems\n', ...
                '----------------------------------------\n', ...
                '\n', ...
                '   (1) 3D Symmetric driven cavity problem. \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
      
      if (~isempty (iopt2))
        switch iopt2
         case 1
          datafile = 'test_stokes_3d_symdrivcav';
          eval (datafile);
        end %# switch (iopt2)
        
        iopt3 = 1; 
        clc;
        fprintf (1, ...
                 ['GeoPDEs_Fluid examples menu: choose function spaces.\n', ...
                  '----------------------------------------------------\n', ...
                  '\n', ...
                  '   (1) TH function spaces.\n \n']);
        
        iopt3 = input ('Please choose a number from above or press <Enter> to return: ');
        
        if (~isempty (iopt3))
          switch iopt3
           case 1
             fprintf (['\nRunning script ' '     ex_bspline_stokes_3d.m', ...
                       '\nwith the data file  %s \n', ...
                       '\nYou can have a look at the source using the command:\n' ...
                       '   type ex_bspline_stokes_3d\n\n'], datafile)

            ex_bspline_stokes_3d;
          end %# switch (iopt3)
          input ('Press <Enter> to continue: ');
        end %# if (~isempty (iopt3))
      end %# while (iopt2 > 0)
    end %# if (~isempty (iopt2))
  end
  
end %# while (iopt > 0)
end
