function DesignSpaceStudy1D(savename,objective_name,server,steps)
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
    Surf = [7, 40];
    Surf_in = [23, 27];
    %savename = 'Obj_WeWm_';
    maxmeshsize_nominal = 0.0015;
    minmeshsize_nominal = 0.00075;
    %objective_name = 'Wm+We';
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename);
    
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
