% GEOPDES_BASE_EXAMPLES: Run some simple examples on how to use the geopdes_base package.
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

function [] = geopdes_base_examples ()

iopt = 1; 
while (iopt > 0)
  clc;
  fprintf (1, ...
           ['GeoPDEs examples menu:\n', ...
            '----------------------\n', ...
            '\n', ...
            '   (1) Examples appearing in the article. \n \n', ...
            '   (2) Other examples in 2D: Poisson problem. \n \n',...
            '   (3) Other examples in 3D: Poisson problem. \n \n']);

  iopt = input ('Please choose a number from above or press <Enter> to return: ');
  clc;


  if (iopt == 1)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs examples menu: examples from the article\n \n', ...
                'C. de Falco, A. Reali, R. Vazquez \n', ...
                'GeoPDEs: a research tool for IsoGeometric Analysis of PDEs (2010) \n', ...
                '------------------------------------------------------------------\n', ...
                '\n', ...
                '   (1) A simple example in 15 lines (section 4).  \n \n', ...
                '   (2) The same example with k-refinement (section 5.1.1). \n \n', ...
                '   (3) The non-isoparametric approach (section 5.1.2). \n \n', ...
                '   (4) A different quadrature rule (section 5.1.3). \n \n', ...
                '   (5) A problem with non-homogeneous boundary conditions (section 5.1.4). \n \n'])
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');
      if (~isempty (iopt2) && iopt2 > 0 && iopt2 < 6)
        vexa = do_example (iopt2);
        disp (vexa);
        eval (vexa);
        input ('Press <Enter> to continue: ');
      end
    end %# while (iopt2 > 0)

  elseif (iopt == 2)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs examples menu: 2D examples\n', ...
                '----------------------------------\n', ...
                '\n', ...
                '   (1) Isoparametric approach: NURBS geometry, NURBS discretization \n \n', ...
                '   (2) Non-isoparametric approach: NURBS geometry, splines discretization \n \n', ...
                '   (3) Laplace eigenvalue problem, non-isoparametric. \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, ...
                   ['GeoPDEs examples menu: 2D examples\n', ...
                    '----------------------------------\n', ...
                    '\n', ...
                    '   2D Poisson problem solved with the isoparametric approach. \n \n', ...
                    '   You can take a look at the script file: examples/ex_nurbs_laplace_2d.m \n \n', ...
                    '   You can solve with one of the following data files: \n \n', ...
                    '   - test_square \n',...
                    '   - test_ring \n',...
                    '   - test_ring_mixed_bc \n',...
                    '   - test_plate_mixed_bc \n \n']);
          fprintf ('Please write the name of the data file: \n');
          data_file = input ('[test_ring_mixed_bc] ', 's');
          if (~isempty (data_file))
            eval(data_file);
          else 
            eval('test_ring_mixed_bc');
          end
          vexa = do_example (6);
          eval (vexa);
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, ...
                   ['GeoPDEs examples menu: 2D examples\n', ...
                    '----------------------------------\n', ...
                    '\n', ...
                    '   2D Poisson problem solved with the non-isoparametric approach. \n \n', ...
                    '   You can take a look at the script file: examples/ex_bspline_laplace_2d.m \n \n', ...
                    '   You can solve with one of the following data files: \n \n', ...
                    '   - test_square \n',...
                    '   - test_ring \n',...
                    '   - test_ring_mixed_bc \n',...
                    '   - test_plate_mixed_bc \n \n']);
          fprintf ('Please write the name of the data file: \n');
          data_file = input ('[test_ring_mixed_bc] ', 's');
          if (~isempty (data_file))
            eval(data_file);
          else 
            eval('test_ring_mixed_bc');
          end
          vexa = do_example (7);
          eval (vexa);
          input ('Press <Enter> to continue: ');

         case 3
          clc;
          fprintf (1, ...
                   ['GeoPDEs examples menu: 2D examples\n', ...
                    '----------------------------------\n', ...
                    '\n', ...
                    '   2D Laplace eigenvalue problem solved with the non-isoparametric approach. \n \n', ...
                    '   You can take a look at the script file: examples/ex_bspline_laplace_2d_eig.m \n \n', ... 
                    '   The only available data file is test_square_eig. \n \n']);
          fprintf ('Please write the name of the data file: \n');
          data_file = input ('[test_square_eig] ', 's');
          if (~isempty (data_file))
            eval(data_file);
          else 
            eval('test_square_eig');
          end
          vexa = do_example (8);
          eval (vexa);
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)
  elseif (iopt == 3)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs examples menu: 3D examples\n', ...
                '----------------------------------\n', ...
                '\n', ...
                '   (1) Isoparametric Poisson problem: NURBS geometry, NURBS discretization \n \n', ...
                '   (2) Non-isoparametric Poisson problem: NURBS geometry, splines discretization \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2))
        switch iopt2
         case 1
          clc;
          fprintf (1, ...
                   ['GeoPDEs examples menu: 3D examples\n', ...
                    '----------------------------------\n', ...
                    '\n', ...
                    '   3D Poisson problem solved with the isoparametric approach. \n \n', ...
                    '   You can take a look at the script file: examples/ex_nurbs_laplace_3d.m \n \n', ...
                    '   You can solve with one of the following data files: \n \n', ...
                    '   - test_cube \n',...
                    '   - test_thick_ring \n \n']);
          fprintf ('Please write the name of the data file: \n');
          data_file = input ('[test_thick_ring] ', 's');
          if (~isempty (data_file))
            eval(data_file);
          else 
            eval('test_thick_ring');
          end
          vexa = do_example (9);
          eval (vexa);
          input ('Press <Enter> to continue: ');
         case 2
          clc;
          fprintf (1, ...
                   ['GeoPDEs examples menu: 3D examples\n', ...
                    '----------------------------------\n', ...
                    '\n', ...
                    '   3D Poisson problem solved with the non-isoparametric approach. \n \n', ...
                    '   You can take a look at the script file: examples/ex_bspline_laplace_3d.m \n \n', ...
                    '   You can solve with one of the following data files: \n \n', ...
                    '   - test_cube \n',...
                    '   - test_thick_ring \n \n']);
          fprintf ('Please write the name of the data file: \n'); 
          data_file = input ('[test_thick_ring] ', 's');
          if (~isempty (data_file))
            eval(data_file);
          else 
            eval('test_thick_ring');
          end
          vexa = do_example (10);
          eval (vexa);
          input ('Press <Enter> to continue: ');
        end %# switch (iopt2)
      end %# if (~isempty (iopt2))
    end %# while (iopt2 > 0)
  end

end %# while (iopt > 0)
end

function vexa = do_example (number)

switch (number)
 case 1
  filename = 'ex_nurbs_laplace_2d_15lines.m';
 case 2
  filename = 'ex_article_section_511.m';
 case 3
  filename = 'ex_article_section_512.m';
 case 4
  filename = 'ex_article_section_513.m';
 case 5
  filename = 'ex_article_section_514.m';
 case 6
  filename = 'ex_nurbs_laplace_2d.m';
 case 7
  filename = 'ex_bspline_laplace_2d.m';
 case 8
  filename = 'ex_bspline_laplace_2d_eig.m';
 case 9
  filename = 'ex_nurbs_laplace_3d.m';
 case 10
  filename = 'ex_bspline_laplace_3d.m';
end

fid = fopen (filename);
vexa = '';
l = fgets (fid);
while (l ~= -1)
  vexa = cat (2, vexa, l);
  l = fgets (fid);
end
end

