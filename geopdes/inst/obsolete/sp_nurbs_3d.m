% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function sp = sp_nurbs_3d (varargin)

warning ('geopdes:obsolete', 'The class SP_NURBS_3D is obsolete. Using SP_NURBS instead')

sp = sp_nurbs (varargin{:});
