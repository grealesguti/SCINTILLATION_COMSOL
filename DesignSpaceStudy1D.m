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
    defaultport = 2036;
    defaultwavelength = '1.5[mm]';
    defaultdeltaY = 0;
    defaultIObject = 'fin';
    defaultWeimag = [];
defaultParts = 2;  % Set default value to 2

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
    addParameter(p, 'port', defaultport, @isnumeric);
    addParameter(p, 'wavelength', defaultwavelength, @ischar);
    addParameter(p, 'deltaY', defaultdeltaY, @isnumeric);
    addParameter(p, 'IObject', defaultIObject, @ischar);
    addParameter(p, 'Weimag', defaultWeimag, @isnumeric);
addParameter(p, 'parts', defaultParts, @isnumeric);  % <-- This line sets 'parts' default to 2
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
    port = p.Results.port;
    wavelength = p.Results.wavelength;
    deltaY = p.Results.deltaY;
    IObject = p.Results.IObject;
    Weimag = p.Results.Weimag;
parts = p.Results.parts;  % <-- Retrieve parts input

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
		mphstart(port);
		pltoption=0;
    else
        fprintf('WARNING: SETTING PLOTS')
		pltoption=1;
	end

    maxmeshsize_nominal = maxmesh;
    minmeshsize_nominal = minmesh;
    
    objectiveFunctionSearch = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surf_in, modelname, pltoption, objective_name, savename, SimpStol,I0,Iend,int,deltaY,IObject,Weimag);
    objectiveFunctionSearch.parts=parts;
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

        figure(3)
        plot(xa(xa > 0), -Rst(xa > 0));
        title(append('Energy Deposition vs Position' , wavelength));
        xlabel('Position (x)');
        ylabel('Energy Deposition (Rst)');
        % Construct full save name
        name_save = fullfile('Rst/', append(savename, objectiveFunctionSearch.creationDate, '.png'));
        % Save the figure
        saveas(gcf, name_save);

    end
% Inside your loop
figure(4)
hold on

    % Plot current curve with label
    plot(xa(xa > 0), -Rst(xa > 0), 'DisplayName', wavelength);
    
    % Update labels and title (safe to call multiple times)
    xlabel('Position (x)');
    ylabel('Energy Deposition (Rst)');
    title('Energy Deposition vs Position for Different Wavelengths');
    
    % Show legend with all plotted curves
    legend('Location', 'best');

% Construct full save name with a unique identifier
name_save = fullfile('Rst/', append(savename, objectiveFunctionSearch.creationDate, '_multiWavelength.png'));

% Save the final figure
saveas(figure(4), name_save);

% --- Figure 5: Normalized curves ---
figure(5);
hold on;

% Extract and normalize data
xvals = xa(xa > 0);
yvals = -Rst(xa > 0);
yvals_norm = (yvals - min(yvals)) / (max(yvals) - min(yvals));  % Normalize to [0, 1]

% Plot normalized curve with wavelength label
plot(xvals, yvals_norm, 'DisplayName',  wavelength);

% Add labels and title
xlabel('Position (x)');
ylabel('Normalized Energy Deposition');
title('Normalized Energy Deposition vs Position for Different Wavelengths');
legend('Location', 'best');

name_save_norm = fullfile('Rst/', append(savename, objectiveFunctionSearch.creationDate, '_normmultiWavelength.png'));

%name_save_norm = fullfile('Rst/', ...
%    sprintf('%s_%s_lambda_%.2fmm_normalized.png', savename, objectiveFunctionSearch.creationDate, '_normmultiWavelength.png'));
saveas(gcf, name_save_norm);

end
