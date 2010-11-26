function pre_install (desc)
  dir =  (fileparts (mfilename ("fullpath")));
  putenv ("GEOPDES_INCLUDE_DIR", [ "-I" canonicalize_file_name(dir) filesep() "inst" filesep() "include"]);
endfunction