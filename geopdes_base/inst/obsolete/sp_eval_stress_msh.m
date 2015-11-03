% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
function [eu, F] = sp_eval_stress_msh (u, space, msh, lambda_lame, mu_lame)

warning ('The function SP_EVAL_STRESS_MSH is obsolete. Using SP_EVAL_MSH instead. Read the help for the usage')

[eu, F] = sp_eval_msh (u, space, msh, 'stress', lambda_lame, mu_lame);