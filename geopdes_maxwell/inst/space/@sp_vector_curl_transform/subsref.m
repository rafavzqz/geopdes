%% override private read access of object data
function value = subsref (obj, S)
value = builtin ('subsref', obj, S);
end