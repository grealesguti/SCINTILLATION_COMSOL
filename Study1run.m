function [Rst,other,model]=Study1run(savename,objective_name,server,x,Ip, varargin)

    fprintf('### 1D OPTIMIZATION MATLAB ###\n')
    % Default values for minmesh and maxmesh
    defaultMinMesh = 0.00075;
    defaultMaxMesh = 0.0015;
    
    if nargin < 5
        fprintf("Number of inputs, %i\n",nargin)
        minmesh = defaultMinMesh;
        fprintf("MinMesh, %f\n",minmesh)
        maxmesh = defaultMaxMesh;
        fprintf("MaxMesh, %f\n",maxmesh)
        SimpStol=1e-6;
        modelname = 'Scintillator3D_1DStudy_2Dgeomv2 - Copy.mph';
    else
        fprintf("Number of inputs, %i\n",nargin)
        minmesh = varargin{1};
        fprintf("MinMesh, %f\n",minmesh)
        maxmesh = varargin{2};
        fprintf("MaxMesh, %f\n",maxmesh)
        SimpStol=varargin{3};
        if strcmp(varargin{4},'spline')
            modelname = 'Scintillator3D_1DStudy_2Dgeom - Spline3.mph';
        else
            modelname = 'Scintillator3D_1DStudy_2Dgeomv2 - Copy.mph';
       end
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
    
    % Define model constants
    if length(x)==1
        Surf = [7, 42];
        Surf_in = [23, 27];
    else
        Surf = [7, 36];
        Surf_in = [23, 19];
    end

    %savename = 'Obj_WeWm_';
    maxmeshsize_nominal = maxmesh;
    minmeshsize_nominal = minmesh;
    %objective_name = 'Wm+We';
    I0=0.005;
    Iend=0.9;
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename, SimpStol,I0,Iend);
    Nvar=1;
    if length(x)==1
        [Rst,other]=objectiveFunctionSearch.RunModel_Case(Ip, x);
    else
        [Rst,other]=objectiveFunctionSearch.RunModel_Case_Nvar(Ip, x, length(x));
    end
    model=objectiveFunctionSearch.model;
  %      Rst = objectiveFunctionSearch.compute(x);
  %      xa = x;
  %      name_save= append('Rst/', savename, objectiveFunctionSearch.creationDate);
  %      saveData(name_save, 'xa', xa, 'Rst', Rst, 'objective_name', objective_name, 'modelname', modelname, 'x0', x0, 'xend', xend, 'steps', steps);
  %      disp(['xa & Rst Variables saved to ', name_save]);
end
