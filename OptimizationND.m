function OptimizationND(savename,objective_name,server, Volcon, optimizer, x0_file)
    % VolCon Default = 1.707000e-04

    LibInitialization()
	if server
		addpath('/apps/generic/comsol/6.2/mli/');
		pause(10);
		mphstart(2036);
		pltoption=0;
	else
		pltoption=1;
	end
    modelname = 'Scintillator3D_1DStudy_2Dgeom - Spline3.mph';
    
    % Define model constants
    Surf = [7, 36];
    Surf_in = [23, 19];
    %savename = 'Obj_WeWm_';
    maxmeshsize_nominal = 0.0015;
    minmeshsize_nominal = 0.00075;
    %objective_name = 'Wm+We';
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename);

    xlim = [0.1,1.9];
    if isempty(x0_file)
        Nvar = 18;
        x0 = ones(Nvar,1)*1; % Default value if x0 file is not provided
    else
        % Read x0 from the provided file
        x0 = load(x0_file);
        Nvar=length(x0);
    end
    
    OBJECTIVE = @(x) objectiveFunctionSearch.compute_Nvar(x,Nvar);
    CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint(x, Volcon, Nvar);
    % Define options for fminsearchbnd
    %if server
    %options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
    %else
    options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
    %end
    varargin = {}; % Define varargin if needed

    %% https://nl.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
    % transformation to a bounded problem using fminsearch. 
    % options and vargaring are directly passed to fminsearch
    % https://nl.mathworks.com/help/matlab/ref/fminsearch.html
    fprintf(append('optimizer:',optimizer))    
    % Select optimizer based on specified case
   switch optimizer
        case 'fminsearchcon'
            % Define options for fminsearchcon
            options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
            varargin = {}; % Define varargin if needed
            % fminsearchcon optimization
            [x_opt, fval, exitflag, output] = fminsearchcon(OBJECTIVE, x0, ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), [], [], CONSTRAINT, options, varargin);
        case 'ga'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            % Genetic algorithm options
            ga_options = optimoptions(@ga, 'Display', 'iter', 'MaxGenerations', 100, 'PopulationSize', 50);
            [x_opt, fval, exitflag, output] = ga(OBJECTIVE, Nvar, [], [], [], [], ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), CONSTRAINT, ga_options);
        case 'patternsearch'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            options = optimoptions(@patternsearch, 'Display', 'iter');
            [x_opt, fval, exitflag, output] = patternsearch(OBJECTIVE, x0, [], [], [], [], ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), CONSTRAINT, options);
        case 'simulannealbnd'
            % NO CONSTRAINT
            options = optimoptions(@simulannealbnd, 'Display', 'iter');
            [x_opt, fval, exitflag, output] = simulannealbnd(OBJECTIVE, x0, ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), options);
        case 'particleswarm'
            % no constraint
            options = optimoptions(@particleswarm, 'Display', 'iter');
            [x_opt, fval, exitflag, output] = particleswarm(OBJECTIVE, Nvar, ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), options);
        case 'fmincon'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            options = optimoptions(@fmincon, 'Display', 'iter');
            [x_opt, fval, exitflag, output] = fmincon(OBJECTIVE, x0, [], [], [], [], ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), CONSTRAINT, options);
        case 'bayesopt'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            % Define the search space
            lb = ones(Nvar, 1) * xlim(1);
            ub = ones(Nvar, 1) * xlim(2);
            search_space = optimizableVariable('x', lb, ub);
            
            % Bayesian optimization options
            bayesopt_options = bayesopt('Display', 'iter');
            [x_opt, fval, exitflag, output] = bayesopt(OBJECTIVE, search_space, 'Constraint', CONSTRAINT, bayesopt_options);
        otherwise
            error('Invalid optimizer specified.');
    end


    name_save= append('Rst/', savename, objectiveFunctionSearch.creationDate);
    saveData(name_save, 'x_opt', x_opt, 'fval', fval, 'exitflag', exitflag, 'output', output);
    disp(['x_opt & fval Variables saved to ', name_save]);
    fprintf('Optimal solution x_opt: %e, fval: %e', x_opt, fval);
end
