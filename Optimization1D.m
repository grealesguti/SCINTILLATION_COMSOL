function Optimization1D(savename, objective_name, server, optimizer, varargin)
    % Create a global structure to store iteration data
    global iterData;
    iterData.iteration = [];
    iterData.funccount = [];
    iterData.fval = [];
    iterData.x = [];



    %DesignSpaceStudy1D('savefile', 'We', 0, 100, 'minmesh', 0.00075, 'maxmesh', 0.0015, 'SimpStol', 1e-5, 'I0', 0.005, 'Iend', 0.9, 'Ampl', 1E7, 'x0',0.75,'xend',1.25);

    fprintf('### 1D OPTIMIZATION MATLAB ###\n');
    
    %Default values for optional parameters
    defaultMinMesh = 0.00075;
    defaultMaxMesh = 0.0015;
    defaultSimpStol = 1e-6;
    defaultI0 = 0.005;
    defaultIend = 0.9;
    defaultAmpl = 1e7;
    defaultx0 = 0.75;
    defaultxend = 1.25;
    defaulttr = 0.75;
    defaulttd = 1.25;
    defaultSurf = [11, 40];
    defaultSurfin =  [23, 27];
    defaultjRINDEX = 8.3559e-06;
    defaultjRINDEX_G = 2.5710e-06;
    defaultjRINDEX_R = 4.1779e-07;
    defaultint = '1D';
    defaultwavelength = '1.5[mm]';
    defaultdeltaY = 0;
    defaultOptx0 = 1;
        defaultIObject = 'fin';


    defaultmodelname = 'Scintillator3D_1DStudy_2Dgeomv2 - Copyv2.mph';
    %Create input parser
    p = inputParser;
    
    %Add required parameters
    addRequired(p, 'savename', @ischar);
    addRequired(p, 'objective_name', @ischar);
    addRequired(p, 'server', @isnumeric);
        addRequired(p, 'optimizer', @ischar);

    %Add optional parameters with default values
    addParameter(p, 'minmesh', defaultMinMesh, @isnumeric);
    addParameter(p, 'maxmesh', defaultMaxMesh, @isnumeric);
    addParameter(p, 'SimpStol', defaultSimpStol, @isnumeric);
    addParameter(p, 'I0', defaultI0, @isnumeric);
    addParameter(p, 'Iend', defaultIend, @isnumeric);
    addParameter(p, 'Ampl', defaultAmpl, @isnumeric);
    addParameter(p, 'x0', defaultx0, @isnumeric);
    addParameter(p, 'xend', defaultxend, @isnumeric);
    addParameter(p, 'tr', defaulttr, @isnumeric);
    addParameter(p, 'td', defaulttd, @isnumeric);
    addParameter(p, 'Surf', defaultSurf, @isnumeric);
    addParameter(p, 'Surfin', defaultSurfin, @isnumeric);
    addParameter(p, 'modelname', defaultmodelname, @ischar);
    addParameter(p, 'jRINDEX', defaultjRINDEX, @isnumeric);
    addParameter(p, 'jRINDEX_G', defaultjRINDEX_G, @isnumeric);
    addParameter(p, 'jRINDEX_R', defaultjRINDEX_R, @isnumeric);
    addParameter(p, 'int', defaultint, @ischar);
    addParameter(p, 'wavelength', defaultwavelength, @ischar);
    addParameter(p, 'deltaY', defaultdeltaY, @isnumeric);
    addParameter(p, 'Optx0', defaultOptx0, @isnumeric);
    addParameter(p, 'IObject', defaultIObject, @ischar);

    %Parse inputs
    parse(p, savename, objective_name, server,optimizer, varargin{:});
    
    %Retrieve values
    minmesh = p.Results.minmesh;
    maxmesh = p.Results.maxmesh;
    SimpStol = p.Results.SimpStol;
    I0 = p.Results.I0;
    Iend = p.Results.Iend;
    Ampl = p.Results.Ampl;
    x0 = p.Results.x0;
    xend = p.Results.xend;    
    tr = p.Results.tr;
    td = p.Results.td;   
    Surf = p.Results.Surf;
    Surf_in = p.Results.Surfin;   
    modelname = p.Results.modelname;   
    jRINDEX = p.Results.jRINDEX;
    jRINDEX_G = p.Results.jRINDEX_G;
    jRINDEX_R = p.Results.jRINDEX_R;
    int = p.Results.int;
    wavelength = p.Results.wavelength;
    deltaY = p.Results.deltaY;
    Optx0 = p.Results.Optx0;
    IObject = p.Results.IObject;

    % Display input values
    fprintf('Number of inputs, %i\n', nargin);
    fprintf('MinMesh, %f\n', minmesh);
    fprintf('MaxMesh, %f\n', maxmesh);
    fprintf('SimpStol, %f\n', SimpStol);
    fprintf('I0, %f\n', I0);
    fprintf('Iend, %f\n', Iend);
    fprintf('Ampl, %f\n', Ampl);
    fprintf('x0, %f\n', x0);
    fprintf('xend, %f\n', xend);
    fprintf('tr, %f\n', tr);
    fprintf('td, %f\n', td);
    fprintf('Surf, [%d, %d]\n', Surf);
    fprintf('Surfin, [%d, %d]\n', Surf_in);
    fprintf('modelname, %s\n', modelname);


    LibInitialization()
	if server
		addpath('/apps/generic/comsol/6.2/mli/');
		pause(10);
		mphstart(2036);
		pltoption=0;
	else
		pltoption=1;
	end
    
    %savename = 'Obj_WeWm_';
    maxmeshsize_nominal = maxmesh;
    minmeshsize_nominal = minmesh;
    %objective_name = 'Wm+We';
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename, SimpStol,I0,Iend,int,deltaY,IObject,[]);
    objectiveFunctionSearch.model.param.set('Ampl', Ampl);
    %objectiveFunctionSearch.model.param.set('t_r', tr);
    %objectiveFunctionSearch.model.param.set('t_d', td);

    % Try setting the new parameters and catch errors
    try
        objectiveFunctionSearch.model.param.set('jRINDEX', jRINDEX);
        objectiveFunctionSearch.model.param.set('jRINDEX_G', jRINDEX_G);
        objectiveFunctionSearch.model.param.set('jRINDEX_R', jRINDEX_R);
        objectiveFunctionSearch.model.param.set('wavelength', wavelength);
    catch ME
        fprintf('Error setting jRINDEX parameters: %s\n', ME.message);
    end
    xlim = [0.1,1.9];
    
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
    filename_optim = append(objectiveFunctionSearch.savename,objectiveFunctionSearch.creationDate,'_iter.mat');
   switch optimizer
        case 'fminsearchcon'
            % Define options for fminsearchcon
            options = optimset('OutputFcn', @(x, optimValues, state) outfun(x, optimValues, state, filename_optim), 'Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
            varargin = {}; % Define varargin if needed
            % fminsearchcon optimization
            [x_opt, fval, exitflag, output] = fminsearchcon(OBJECTIVE, Optx0, ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), [], [], [], options, varargin);
        case 'ga'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            % Genetic algorithm options
            ga_options = optimoptions(@ga, 'Display', 'iter', 'MaxGenerations', 100, 'PopulationSize', 50, 'OutputFcn', @(x, optimValues, state) outfun(x, optimValues, state, filename_optim));
            [x_opt, fval, exitflag, output] = ga(OBJECTIVE, Nvar, [], [], [], [], ones(Nvar,1)*xlim(1), ones(Nvar,1)*xlim(2), [], ga_options);
        case 'patternsearch'
            CONSTRAINT = @(x) objectiveFunctionSearch.ComsolVolumeConstraint_ga(x, Volcon, Nvar);
            options = optimoptions(@patternsearch, 'Display', 'iter','OutputFcn', @(x, optimValues, state) outfun(x, optimValues, state, filename_optim));
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
            options = optimoptions(@fmincon, 'Display', 'iter',...
                'OutputFcn', @(x, optimValues, state) outfun(x, optimValues, state, filename_optim),...
                'Algorithm', 'sqp',...
                'HessianApproximation', 'bfgs');
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
