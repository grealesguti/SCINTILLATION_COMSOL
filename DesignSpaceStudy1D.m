function DesignSpaceStudy1D(savename,objective_name,server,steps, varargin)


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
    defaultSurfin =  [21, 25];
    defaultmodelname = 'Scintillator3D_1DStudy_2Dgeomv2 - Copyv2.mph';
    %Create input parser
    p = inputParser;
    
    %Add required parameters
    addRequired(p, 'savename', @ischar);
    addRequired(p, 'server', @isnumeric);
    addRequired(p, 'steps', @isnumeric);
    
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
    addParameter(p, 'modelname', defaultmodelname, @isnumeric);

    %Parse inputs
    parse(p, savename, server, steps, varargin{:});
    
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

    %Display input values
    fprintf('Number of inputs, %i\n', nargin);
    fprintf('MinMesh, %f\n', minmesh);
    fprintf('MaxMesh, %f\n', maxmesh);
    fprintf('SimpStol, %f\n', SimpStol);
    fprintf('I0, %f\n', I0);
    fprintf('Iend, %f\n', Iend);

    LibInitialization()
	if server
		addpath('/apps/generic/comsol/6.2/mli/');
		pause(10);
		mphstart(2036);
		pltoption=0;
	else
		pltoption=1;
	end

    maxmeshsize_nominal = maxmesh;
    minmeshsize_nominal = minmesh;
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename, SimpStol,I0,Iend);
    objectiveFunctionSearch.model.param.set('Ampl', Ampl);


    steps = 50;
    Rst = zeros(steps, 1);
    xa = zeros(steps, 1);
    
    for i = 1:steps
        x = x0 + (xend - x0) / (steps-1) * (i - 1);
        disp(x)
        Rst(i) = objectiveFunctionSearch.compute(x);
        xa(i) = x;
        name_save= append('Rst/', savename, objectiveFunctionSearch.creationDate);
        saveData(name_save, 'xa', xa, 'Rst', Rst, 'objective_name', objective_name, 'modelname', modelname, 'x0', x0, 'xend', xend, 'steps', steps);
        disp(['xa & Rst Variables saved to ', name_save]);
    end
end
