% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [eu, F] = sp_eval_div_msh (u, space, msh)

warning ('geopdes:obsolete', 'The function SP_EVAL_DIV_MSH is obsolete. Using SP_EVAL_MSH instead. Read the help for the usage')

[eu, F] = sp_eval_msh (u, space, msh, 'divergence');
