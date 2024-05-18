function stop = outfun(x, optimValues, state, filename)
    global iterData;

    stop = false; % Do not stop the algorithm

    switch state
        case 'init'
            % Initialization state
            iterData.iteration = [];
            iterData.funccount = [];
            iterData.fval = [];
            iterData.x = [];
        case 'iter'
            % Collect data at each iteration
            if isfield(optimValues, 'iteration') % for fmincon and fminunc
                iter = optimValues.iteration;
            elseif isfield(optimValues, 'funccount') % for patternsearch and ga
                iter = optimValues.funccount;
            else
                iter = NaN; % fallback if neither field is found
            end
            
            iterData.iteration = [iterData.iteration; iter];
            iterData.funccount = [iterData.funccount; optimValues.funccount];
            iterData.fval = [iterData.fval; optimValues.fval];
            iterData.x = [iterData.x; x'];
            % Save the iteration data to the specified .mat file
            save(filename, 'iterData');
            % Print the latest iteration data
            fprintf('Iteration: %d, Function Count: %d, Function Value: %.6f, x: %s\n', ...
                    iter, optimValues.funccount, optimValues.fval, x);
        case 'done'
            % Final state
            % Optionally save the final iteration data
            save(filename, 'iterData');
    end
end
