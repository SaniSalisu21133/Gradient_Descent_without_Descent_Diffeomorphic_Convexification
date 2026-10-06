function k = style_key(name)
%STYLE_KEY  Map a run name to a field of STYLE_DEFS (field names cannot contain '-').
k = strrep(name, '-', '');
end
