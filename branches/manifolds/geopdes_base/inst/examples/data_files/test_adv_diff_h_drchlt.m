function gh = test_adv_diff_h_drchlt (x, y, ind)

%index refers to the geopdes numbering of the sides
  switch (ind)
    case 1
      gh = zeros (size(x));
    case 2
      gh = y>0.8;
    case 3
      gh = zeros (size(x));
    case 4
      gh = ones (size(x));
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end

