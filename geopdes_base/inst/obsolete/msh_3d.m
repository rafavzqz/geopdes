% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function msh = msh_3d (breaks, qn, qw, geo, varargin)

warning ('geopdes:obsolete', 'The class MSH_3D is obsolete. Using MSH_CARTESIAN instead')

msh = msh_cartesian (breaks, qn, qw, geo, varargin{:});
