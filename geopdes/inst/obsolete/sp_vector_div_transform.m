% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_vector_div_transform (scalar_spaces, msh)

  warning ('geopdes:obsolete', 'The function SP_VECTOR_DIV_TRANSFORM is obsolete. Using SP_VECTOR with the transformation as last argument')
  sp = sp_vector (scalar_spaces, msh, 'div-preserving');

end