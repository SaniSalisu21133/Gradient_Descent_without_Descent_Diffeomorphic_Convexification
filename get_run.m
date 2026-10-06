function r = get_run(runs, name)
%GET_RUN  Return the run with the given name from a cell array of run records.
for k = 1:numel(runs)
    if strcmp(runs{k}.name, name)
        r = runs{k};  return
    end
end
error('run %s not found', name);
end
