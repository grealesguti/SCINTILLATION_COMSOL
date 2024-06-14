function DesignSpaceStudy1D(savename,objective_name,server,steps, varargin)


    %DesignSpaceStudy1D('savefile', 'We', 0, 100, 'minmesh', 0.00075, 'maxmesh', 0.0015, 'SimpStol', 1e-5, 'I0', 0.005, 'Iend', 0.9, 'Ampl', 1E7, 'x0',0.75,'xend',1.25);

    fprintf('### Study1D MATLAB ###\n');
    
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

    defaultmodelname = 'Scintillator3D_1DStudy_2Dgeomv2 - Copyv2.mph';
    %Create input parser
    p = inputParser;
    
    %Add required parameters
    addRequired(p, 'savename', @ischar);
    addRequired(p, 'objective_name', @ischar);
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
    addParameter(p, 'modelname', defaultmodelname, @ischar);
    addParameter(p, 'jRINDEX', defaultjRINDEX, @isnumeric);
    addParameter(p, 'jRINDEX_G', defaultjRINDEX_G, @isnumeric);
    addParameter(p, 'jRINDEX_R', defaultjRINDEX_R, @isnumeric);
    addParameter(p, 'int', defaultint, @ischar);

    %Parse inputs
    parse(p, savename,objective_name, server, steps, varargin{:});
    
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


    %Display input values
    fprintf('Number of inputs, %i\n', nargin);
    fprintf('MinMesh, %f\n', minmesh);
    fprintf('MaxMesh, %f\n', maxmesh);
    fprintf('SimpStol, %f\n', SimpStol);
    fprintf('I0, %f\n', I0);
    fprintf('Iend, %f\n', Iend);

    LibInitialization()
	if server
        fprintf('Starting Server')
		addpath('/apps/generic/comsol/6.2/mli/');
		pause(10);
		mphstart(2036);
		pltoption=0;
    else
        fprintf('WARNING: SETTING PLOTS')
		pltoption=1;
	end

    maxmeshsize_nominal = maxmesh;
    minmeshsize_nominal = minmesh;
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename, SimpStol,I0,Iend,int);
    objectiveFunctionSearch.model.param.set('Ampl', Ampl);
    %objectiveFunctionSearch.model.param.set('t_r', tr);
    %objectiveFunctionSearch.model.param.set('t_d', td);

    % Try setting the new parameters and catch errors
    try
        objectiveFunctionSearch.model.param.set('jRINDEX', jRINDEX);
        objectiveFunctionSearch.model.param.set('jRINDEX_G', jRINDEX_G);
        objectiveFunctionSearch.model.param.set('jRINDEX_R', jRINDEX_R);
    catch ME
        fprintf('Error setting jRINDEX parameters: %s\n', ME.message);
    end

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
