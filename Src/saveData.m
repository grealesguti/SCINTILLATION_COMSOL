function saveData(filename, varargin)
    % Function to save variables to a .mat file
    % Usage: saveVariablesToMatFile(filename, var1, var2, ...)
    %   filename: Name of the .mat file to save
    %   varargin: List of variables to be saved
    
    % Create a structure to store the variables
    data = struct();
    
    % Assign each input variable to the structure
    for i = 1:numel(varargin)/2
        varname = varargin{((i-1)*2)+1}; % Get variable name
        data.(varname) = varargin{((i-1)*2)+2}; % Assign variable to structure
    end
    
    % Save to .mat file
    save(fullfile([filename, '.mat']), '-struct', 'data');
    
    % Save to .csv files
    fields = fieldnames(data);
    for i = 1:numel(fields)
        csvwrite(fullfile([filename, '_', fields{i}, '.csv']), data.(fields{i}));
    end


end
