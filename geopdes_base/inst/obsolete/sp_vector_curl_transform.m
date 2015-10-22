% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_vector_curl_transform (scalar_spaces, msh)

  warning ('geopdes:obsolete', 'The function SP_VECTOR_CURL_TRANSFORM is obsolete. Using SP_VECTOR with the transformation as last argument')
  sp = sp_vector (scalar_spaces, msh, 'curl-preserving');

end