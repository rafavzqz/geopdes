% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function varargout = geo_2d_nurbs (nurbs, pts, ders)

warning ('geopdes:obsolete', 'The function GEO_2D_NURBS is obsolete. Using GEO_NURBS instead')

varargout = geo_nurbs (nurbs, pts, ders);
