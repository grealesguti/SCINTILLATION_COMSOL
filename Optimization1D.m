function Optimization1D(savename,objective_name,server,optimizer,varargin)

    % Default values for minmesh and maxmesh
    defaultMinMesh = 0.00075;
    defaultMaxMesh = 0.0015;
    
    % Check if minmesh and maxmesh are provided, otherwise use defaults
    if nargin < 5
        fprintf("Number of inputs, %i",nargin)
        minmesh = defaultMinMesh;
        fprintf("MinMesh, %f",minmesh)
        maxmesh = defaultMaxMesh;
        fprintf("MaxMesh, %f",maxmesh)
    else
        fprintf("Number of inputs, %i",nargin)
        minmesh = varargin{1};
        fprintf("MinMesh, %f",minmesh)
        maxmesh = varargin{2};
        fprintf("MaxMesh, %f",maxmesh)
    end
    fprintf(append('Objective name: ',objective_name))
    fprintf(append('Optimizer name: ',optimizer))


    LibInitialization()
	if server
		addpath('/apps/generic/comsol/6.2/mli/');
		pause(10);
		mphstart(2036);
		pltoption=0;
	else
		pltoption=1;
	end
    modelname = 'Scintillator3D_1DStudy_2Dgeomv2 - Copy.mph';
    
    % Define model constants
    Surf = [7, 42];
    Surf_in = [23, 27];
    %savename = 'Obj_WeWm_';
    maxmeshsize_nominal = maxmesh;
    minmeshsize_nominal = minmesh;
    %objective_name = 'Wm+We';
    
   objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename);

    xlim = [0.1,1.9];
    x0=1;
    
    OBJECTIVE = @(x) objectiveFunctionSearch.compute(x);
    
    % Define options for fminsearchbnd
    if server
    options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
    else
    options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
    end
    varargin = {}; % Define varargin if needed

    % https://nl.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
    % transformation to a bounded problem using fminsearch. 
    % options and vargaring are directly passed to fminsearch
    % https://nl.mathworks.com/help/matlab/ref/fminsearch.html
   % [x_opt,fval,exitflag,output] = fminsearchbnd(OBJECTIVE,x0,xlim(1),xlim(2),options,varargin);

    fprintf(append('optimizer:',optimizer))    
    % Select optimizer based on specified case
    Nvar=1;
   switch optimizer
        case 'fminsearchcon'
            % Define options for fminsearchcon
            options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
            varargin = {}; % Define varargin if needed
            % fminsearchcon optimization
            [x_opt, fval, exitflag, output] = fminsearchcon(OBJECTIVE, x0, ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), [], [], [], options, varargin);
        case 'ga'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            % Genetic algorithm options
            ga_options = optimoptions(@ga, 'Display', 'iter', 'MaxGenerations', 100, 'PopulationSize', 50);
            [x_opt, fval, exitflag, output] = ga(OBJECTIVE, Nvar, [], [], [], [], ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), [], ga_options);
        case 'patternsearch'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            options = optimoptions(@patternsearch, 'Display', 'iter');
            [x_opt, fval, exitflag, output] = patternsearch(OBJECTIVE, x0, [], [], [], [], ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), [], options);
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
            [x_opt, fval, exitflag, output] = fmincon(OBJECTIVE, x0, [], [], [], [], ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), [], options);
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
