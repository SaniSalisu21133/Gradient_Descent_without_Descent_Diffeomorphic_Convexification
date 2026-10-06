function save_fig(fig, name)
%SAVE_FIG  Export a figure as vector PDF and as 300-dpi PNG (MATLAB R2020a or later).
%   The warning about long vectorization times for dense contour plots is harmless and is
%   switched off here; the PDF stays fully vectorial.
w = warning('off', 'all');
drawnow;
exportgraphics(fig, [name '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
exportgraphics(fig, [name '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
warning(w);
fprintf('saved %s.pdf / .png\n', name);
end
