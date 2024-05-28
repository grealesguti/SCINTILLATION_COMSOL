function DesignSpaceStudy1D(savename,objective_name,server,steps, varargin)

    fprintf('### 1D OPTIMIZATION MATLAB ###\n')
    % Default values for minmesh and maxmesh
    defaultMinMesh = 0.00075;
    defaultMaxMesh = 0.0015;
    defaultI0=0.005;
    defaultIend=0.9;
    
    if nargin < 5
        fprintf("Number of inputs, %i\n",nargin)
        minmesh = defaultMinMesh;
        fprintf("MinMesh, %f\n",minmesh)
        maxmesh = defaultMaxMesh;
        fprintf("MaxMesh, %f\n",maxmesh)
        SimpStol=1e-6;
        I0=defaultI0;
        Iend=defaultIend;
    else
        fprintf("Number of inputs, %i\n",nargin)
        minmesh = varargin{1};
        fprintf("MinMesh, %f\n",minmesh)
        maxmesh = varargin{2};
        fprintf("MaxMesh, %f\n",maxmesh)
        SimpStol=varargin{3};
        I0=varargin{4};
        Iend=varargin{5};

    end
    fprintf(append('Objective name: ',objective_name,'\n'))

    LibInitialization()
	if server
		addpath('/apps/generic/comsol/6.2/mli/');
		pause(10);
		mphstart(2036);
		pltoption=0;
	else
		pltoption=1;
	end
    modelname = 'Scintillator3D_1DStudy_2Dgeomv2 - Copyv2.mph';
    
    % Define model constants
    Surf = [11, 38];
    Surf_in = [21, 25];
    %savename = 'Obj_WeWm_';
    maxmeshsize_nominal = maxmesh;
    minmeshsize_nominal = minmesh;
    %objective_name = 'Wm+We';
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename, SimpStol,I0,Iend);
    
    x0 = 0.1;
    xend = 1.9;
    %steps = 50;
    Rst = zeros(steps + 2, 1);
    xa = zeros(steps + 2, 1);
    
    for i = 1:steps
        x = x0 + (xend - x0) / steps * (i - 1);
        disp(x)
        Rst(i) = objectiveFunctionSearch.compute(x);
        xa(i) = x;
        name_save= append('Rst/', savename, objectiveFunctionSearch.creationDate);
        saveData(name_save, 'xa', xa, 'Rst', Rst, 'objective_name', objective_name, 'modelname', modelname, 'x0', x0, 'xend', xend, 'steps', steps);
        disp(['xa & Rst Variables saved to ', name_save]);
    end
end
