% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function msh = msh_2d (breaks, qn, qw, geo, varargin)

warning ('The class MSH_2D is obsolete. Using MSH_CARTESIAN instead')

msh = msh_cartesian (breaks, qn, qw, geo, varargin{:});
