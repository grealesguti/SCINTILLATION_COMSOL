function Optimization1D(savename,objective_name,server,varargin)

    % Default values for minmesh and maxmesh
    defaultMinMesh = 0.00075;
    defaultMaxMesh = 0.0015;
    
    % Check if minmesh and maxmesh are provided, otherwise use defaults
    if nargin < 4
        minmesh = defaultMinMesh;
        maxmesh = defaultMaxMesh;
    else
        minmesh = varargin{1};
        maxmesh = varargin{2};
    end

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
    [x_opt,fval,exitflag,output] = fminsearchbnd(OBJECTIVE,x0,xlim(1),xlim(2),options,varargin);

    name_save= append('Rst/', savename, objectiveFunctionSearch.creationDate);
    saveData(name_save, 'x_opt', x_opt, 'fval', fval, 'exitflag', exitflag, 'output', output);
    disp(['x_opt & fval Variables saved to ', name_save]);
    fprintf('Optimal solution x_opt: %e, fval: %e', x_opt, fval);
end
