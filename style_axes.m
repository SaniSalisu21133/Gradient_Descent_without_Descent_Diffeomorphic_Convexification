function style_axes(ax, ttl)
%STYLE_AXES  Boxed axes, light dotted grid and LaTeX labels.
S = style_defs();
set(ax, 'FontSize', 10, 'Box', 'on', 'TickDir', 'in', 'LineWidth', 0.9, 'Layer', 'top', ...
        'TickLabelInterpreter', 'latex', 'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', ':', 'GridColor', S.grid, 'GridAlpha', 1, ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
if nargin > 1 && ~isempty(ttl)
    title(ax, ttl, 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 11);
end
end
