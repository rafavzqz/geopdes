% msh_partition: partititon a mesh in parts suitable for parallel execution
%
%  [breaks boundaries] = msh_partition( pieces, data ...)
%
% INPUT:
%
% pieces:   number of required parts
% data:	    a list of data in the [data, value] format. Possible data are
%               'breaks', partition trying to divide the elements evenly
%               'knots', if also 'degree' is present partition trying to distribute 
%                        degrees of freedom evenly, else fallback to distributing elements
%               'degree', gives the degree of the spline functions
%               'index', if set return onl the breaks for the part having that label else return a
%                        cell array with all the parts.
%
% OUTPUT:
%
% breaks:     a cell array containing the breaks along all dimension of all the parts
% boundaries: a cell array containing for each part the logical vector of
%             needed boundaries
%
% Copyright (C) 2010 Andrea Bressan
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [breaks boundaries] = msh_partition( pieces, varargin)

has_breaks = false;
has_knots  = false;
has_degree = false;
has_index  = false;

if (rem(length(varargin),2))
  error('msh_partition: data must be given in the [name, value] format');
else
  for ii=1:2:length(varargin)
    if strcmp(varargin{ii}, 'breaks')
      full_breaks=varargin{ii+1};
      has_breaks=true;
    elseif strcmp(varargin{ii}, 'knots')
      knots=varargin{ii+1};
      has_knots=true;
    elseif strcmp(varargin{ii}, 'degree')
      degree=varargin{ii+1}; 
      has_degree=true;
    elseif strcmp(varargin{ii}, 'index')
      index=varargin{ii+1};
      has_index=true;
    end
  end
end

if has_index && (index<=0 || index>pieces)
   error('msh_partition: invalid value for parameter ''index'', should be between 1 and %d', pieces);
end 
if has_breaks && has_knots
  warning('msh_partition: both ''breaks'' and ''knots'' given, knots are used.');
end
if ~(has_breaks || has_knots)
  error('msh_partition: at least ''breaks'' or ''knots'' should be provided');
end
if has_knots
  for ii=1:length(knots)
    full_breaks{ii}=unique(knots{ii});
  end
end
if (has_knots && has_degree)
  if (length(knots)~= length(degree))
    error('msh_partition: knots and degree lengths disagree');
  else
    [breaks boundaries] = partition_dofs (pieces, full_breaks,knots,degree);
  end
else
  [breaks boundaries] = partition_elements (pieces, full_breaks);
end
if length(breaks)<pieces 
  if ~ has_index
    warning('msh_partition: cannot partition in %d pieces; only %d used',pieces, prod(divisions));
  else
    error('msh_partition: cannot partition in %d pieces; only %d used',pieces, prod(divisions));
  end
end
if (has_index)
  breaks     = breaks{index};
  boundaries = boundaries{index};
end
end

function [breaks boundaries] = partition_elements (pieces, full_breaks)
  for i=1:length(full_breaks)
    dim(i)=length(full_breaks{i})-1;
  end
  max_div=dim;
  [cut_points divisions] = partition_equally(pieces,dim,max_div);
  for i=1:prod(divisions)
    part_coord = get_coordinates(divisisions,i);
    for j=1:length(full_breaks)
      part_breaks{j}=full_breaks{j}(cut_points{j}(part_coord(j)):cut_points{j}(part_coord(j)+1));
    end
    breaks{i}  = part_breaks;
    bnd_vector = zeros([1 2*length(dim)]);
    bnd_vector(1:2:2*length(dim)) = part_coord == 1;
    bnd_vector(2:2:2*length(dim)) = part_coord == divisions;
    boundaries{i} = find(bnd_vector);
  end
end
  
function [breaks boundaries] = partition_dofs (pieces, full_breaks,knots,degree)
  for i=1:length(knots)
    dim(i)=length(knots{i})-1-degree(i);
    max_div(i)=length(full_breaks{i})-1;
  end
  max_div=min(max_div,floor(dim./(degree+1)));
  [cut_points divisions] = partition_equally(pieces,dim,max_div);
  for i=1:prod(divisions)
      part_coord = get_coordinates(divisions,i);
    for j=1:length(full_breaks)
      part_start = knots{j}(cut_points{j}(part_coord(j)));
      part_end   = knots{j}(cut_points{j}(part_coord(j)+1));
      part_break_start = find(part_start==full_breaks{j});
      part_break_end   = find(part_end  ==full_breaks{j});
      part_breaks{j}=full_breaks{j}(part_break_start:part_break_end);
    end
    breaks{i}  = part_breaks;
    bnd_vector = zeros([1 2*length(dim)]);
    bnd_vector(1:2:2*length(dim)) = part_coord == 1;
    bnd_vector(2:2:2*length(dim)) = part_coord == divisions;
    boundaries{i} = find(bnd_vector);
  end
end


function [cut_points div] = partition_equally(pieces, dim, max_div)
% a very simple algorithm: factorize the number of piece and divide
% at each step the biggest edge by the biggest factor. If this is 
% not possible try to best approximate the factor.

steps   = factor(pieces);
div     = ones(size(dim)); % div(i) = number of pieces along the i-Th coordinate
bdim    = dim;             % bdim(i) = dimension of a piece along the i-Th
                           %           coordinate
min_dim = dim ./ max_div;  % minimum allowed dimension.

for i=length(steps):-1:1
  [sbdim reordered]=sort(bdim,'descend'); % find biggest edge 
  can_divide = div * steps(i) <= max_div;
  can_divide = can_divide(reordered);
  indices = find(can_divide);
  if ~isempty(indices) % there is a direction in which it is possible to divide
                       % so we do this
    index   = reordered(indices(1));
    div(index) = div(index)*steps(i);
  else    % the factor is too big thus we cannot divide
          % we just try to get near the requested divisions by maximixing divisions
    requested  = prod(steps(i:end));
    [gain pos] = sort(bdim ./ min_dim  ,'descend'); 
    for j=1:length(gain)
      total_gain(j)=prod(gain(1:j));
    end
    indices = find ( total_gain > requested/prod(div));
    if ~isempty(indices)
      last = indices(1)-1;
    else
      last = length(dim);
    end
    for j=1:last
      div(pos(j)) = max_div(pos(j));
    end
  end
  bdim = dim ./ div;    
end
bdim=floor(bdim); % basic part size
reminder = dim - bdim .* div;
regain   = ceil(div ./ reminder);
for i=1:length(dim)
  if regain(i) ~= Inf
    bonus = repmat (1:ceil(reminder(i)+1), [regain(i) 1]);
    bonus = bonus(:);
    bonus = bonus(floor(regain(i)/2):floor(regain(i)/2)+div(i))';
    bonus(div(i)+1) = dim(i)+1 - bdim(i)*div(i);
  else
    bonus = ones([1 div(i)+1]);
  end
  cut_points{i} = (0:div(i))*bdim(i)+bonus;
end
end

function coord = get_coordinates (dims, ind )
  nd = length (dims);
  coord = ones(length(ind),nd);
  scale = [1; cumprod(dims(:))];
  for i = nd:-1:2
    k = (ind >= scale(i));
    r = ones (size (ind));
    t = zeros (size (ind));
    t(k) = floor ((ind(k) - 1) / scale(i));
    r(k) = t(k) + 1;
    coord(:,i) = r;
    ind(k) -= t(k) * scale(i);
  end
  coord(:,1) = ind;
end