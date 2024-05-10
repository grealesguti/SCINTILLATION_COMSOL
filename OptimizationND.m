function OptimizationND(savename,objective_name,server, Volcon)
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
    Surf = [7, 40];
    Surf_in = [23, 27];
    %savename = 'Obj_WeWm_';
    maxmeshsize_nominal = 0.0015;
    minmeshsize_nominal = 0.00075;
    %objective_name = 'Wm+We';
    Nvar = 18;
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename);

    xlim = [0.1,1.9];
    x0=1;
    
    OBJECTIVE = @(x) objectiveFunctionSearch.compute_Nvar(x,Nvar);
    CONSTRAINT = 
    % Define options for fminsearchbnd
    if server
    options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
    else
    options = optimset('OutputFcn', @objectiveFunctionSearch.outfun, 'Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 200, 'MaxFunEvals', 500);
    end
    varargin = {}; % Define varargin if needed

    %% https://nl.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
    % transformation to a bounded problem using fminsearch. 
    % options and vargaring are directly passed to fminsearch
    % https://nl.mathworks.com/help/matlab/ref/fminsearch.html
    
    % MISSING VOLUME CONSTRAINT!!!
    [x_opt,fval,exitflag,output] = fminsearchcon(OBJECTIVE,ones(Nvar,1)*x0,ones(Nvar,1)*xlim(1),ones(Nvar,1)*xlim(2),options,varargin);

    name_save= append('Rst/', savename, objectiveFunctionSearch.creationDate);
    saveData(name_save, 'x_opt', x_opt, 'fval', fval, 'exitflag', exitflag, 'output', output);
    disp(['x_opt & fval Variables saved to ', name_save]);
    fprintf('Optimal solution x_opt: %e, fval: %e', x_opt, fval);
end
