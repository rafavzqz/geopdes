% GEOPDES_BASE_EXAMPLES: Run some simple examples on how to use the geopdes_base package.
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
                '  2D Poisson problem solved with the isoparametric approach: \n' ...
                '   NURBS geometry, NURBS discretization \n \n', ...
                '       (1) Plate with a hole. \n', ...
                '       (2) Quarter of a ring with Dirichlet boundary conditions. \n', ...
                '       (3) Quarter of a ring with mixed boundary conditions. \n \n', ...
                '  2D Poisson problem solved with the non-isoparametric approach: \n' ...
                '   NURBS (or splines) geometry, B-splines discretization \n \n', ...
                '       (4) Unit square domain. \n', ...
                '       (5) Plate with a hole. \n', ...
                '       (6) Quarter of a ring with Dirichlet boundary conditions. \n', ...
                '       (7) Quarter of a ring with mixed boundary conditions. \n \n', ...
                '   2D Laplace eigenvalue problem, non-isoparametric \n \n', ...
                '       (8) Unit square domain. \n \n']);

      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2) && iopt2 > 0 && iopt2 < 9)
        [vexa, filename] = do_example (iopt2+5);
        clc
        fprintf (1, 'You can have a look at the source file: %s \n \n', filename);
        eval (vexa);
        input ('Press <Enter> to continue: ');
      end
    end  %# while (iopt2 > 0)

  elseif (iopt == 3)
    iopt2 = 1; 
    while (iopt2 > 0)
      clc;
      fprintf (1, ...
               ['GeoPDEs examples menu: 3D examples\n', ...
                '----------------------------------\n', ...
                '\n', ...
                '  3D Poisson problem solved with the isoparametric approach: \n', ...
                '   NURBS geometry, NURBS discretization \n \n', ...
                '       (1) Thick quarter of a ring. \n \n', ...
                '  3D Poisson problem solved with the non-isoparametric approach: \n', ...
                '   NURBS (or splines) geometry, B-splines discretization \n \n', ...
                '       (2) Unit cube domain. \n \n', ...
                '       (3) Thick quarter of a ring. \n \n']);
      
      iopt2 = input ('Please choose a number from above or press <Enter> to return: ');

      if (~isempty (iopt2) && iopt2 > 0 && iopt2 < 4)
        [vexa, filename] = do_example (13+iopt2);
        clc
        fprintf (1, 'You can have a look at the source file: %s \n \n', filename);        
        eval (vexa);
        input ('Press <Enter> to continue: ');
      end
    end %# while (iopt2 > 0)
  end

end %# while (iopt > 0)
end

function [vexa, filename] = do_example (number)

switch (number)
 case 1
  filename = 'ex_article_15lines.m';
 case 2
  filename = 'ex_article_section_511.m';
 case 3
  filename = 'ex_article_section_512.m';
 case 4
  filename = 'ex_article_section_513.m';
 case 5
  filename = 'ex_article_section_514.m';
 case 6
  filename = 'ex_laplace_iso_plate.m';
 case 7
  filename = 'ex_laplace_iso_ring.m';
 case 8
  filename = 'ex_laplace_iso_ring_mixed_bc.m';
 case 9
  filename = 'ex_laplace_square.m';
 case 10
  filename = 'ex_laplace_plate.m';
 case 11
  filename = 'ex_laplace_ring.m';
 case 12
  filename = 'ex_laplace_ring_mixed_bc.m';
 case 13
  filename = 'ex_laplace_eig_square.m';
 case 14
  filename = 'ex_laplace_iso_thick_ring.m';  
 case 15
  filename = 'ex_laplace_cube.m';  
 case 16
  filename = 'ex_laplace_thick_ring.m';  
end

fid = fopen (filename);
vexa = '';
l = fgets (fid);
while (l ~= -1)
  vexa = cat (2, vexa, l);
  l = fgets (fid);
end
end

