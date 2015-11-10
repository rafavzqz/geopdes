%% override private write access of objects in this class
function obj = subsasgn (obj, S, value)
obj = builtin('subsasgn', obj, S, value);
end