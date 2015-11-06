% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_nurbs_2d (varargin)

warning ('geopdes:obsolete', 'The class SP_NURBS_2D is obsolete. Using SP_NURBS instead')

sp = sp_nurbs (varargin{:});
