function o = parse_opts(args, defaults)
%PARSE_OPTS  Minimal name/value parser:  o = parse_opts(varargin, defaults).
o = defaults;
if mod(numel(args), 2) ~= 0
    error('options must be given as name/value pairs');
end
for i = 1:2:numel(args)
    key = args{i};
    if ~isfield(defaults, key)
        error('unknown option "%s"', key);
    end
    o.(key) = args{i+1};
end
end
