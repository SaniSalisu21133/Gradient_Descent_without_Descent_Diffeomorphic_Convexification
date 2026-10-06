function h = add_series(ax, key, x, y, varargin)
%ADD_SERIES  Plot one method's curve with its fixed colour and line style (no markers).
S = style_defs();
s = S.(style_key(key));
x = x(:);  y = y(:);
ok = isfinite(x) & isfinite(y);
h = plot(ax, x(ok), y(ok), 'Color', s.color, 'LineStyle', s.ls, 'LineWidth', s.lw, ...
         'DisplayName', s.label, varargin{:});
end
