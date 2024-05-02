classdef Objective1DAdapSearch
    properties
        minmeshsize_nominal
        maxmeshsize_nominal
        Surf
        model
        creationDate
        plt
        objective
        Surfin
    end
    
    methods
        function obj = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surfin, model, plt, objective)
            import com.comsol.model.util.* 
            obj.plt = plt;
            if plt
                ModelUtil.showProgress(true);
            end
            obj.minmeshsize_nominal = minmeshsize_nominal;
            obj.maxmeshsize_nominal = maxmeshsize_nominal;
            obj.Surf = Surf;
            obj.Surfin = Surfin;        
            obj.model = mphopen(model);
            obj.creationDate = datestr(now, 'yyyymmdd');
            obj.objective = objective;
        end
        
        function WI = compute(obj, x)
            % Print vector x to the screen
            fprintf('Vector x: ');
            disp(x);

            LYSO_L = 0.0275;
            I0 = 0.01;
            IEND = 0.9;
            x0 = 0.05;
            xend = 1.95;

            if (x0 < x) && (x < xend)
                FUN = @(Ip) obj.ObjectiveQuad_1D(Ip, x);
                [WI, ~, ~] = adaptiveSimpson(FUN, I0, IEND, 'parts', 2);
                WI = -WI;
                disp(['Objective value: ', sprintf('%.2e', WI)]);
            else
                WI = 1000;
            end
        end
        
        function WI = ObjectiveQuad_1D(obj, Ip, x)
            % Initialize Run function
            FUN = @(xx) obj.RunModel_Case(xx, x);
            WI = zeros(length(Ip), 1);
            for i = 1:length(Ip)
                disp(['Solving impact number: ', num2str(i), ' Out of: ', num2str(length(Ip)), '. Location: ', num2str(Ip(i))]);
                WI(i) = FUN(Ip(i));
            end
            WI = WI';
        end
        
        function [Wet,Wmt,We_int,Wm_int,We,Wm,We_in,Wm_in] = RunModel_AllRst(obj, Ip, x)
            mesh_min = obj.minmeshsize_nominal;
            mesh_max = obj.maxmeshsize_nominal;
            obj.model = obj.model;
            obj.Surf = obj.Surf;
            
            obj.model.param.set('maxmeshsize', mesh_max);
            obj.model.param.set('minmeshsize', mesh_min);
            obj.model.param.set('impact_location', Ip); % loaded as percentage over 1
            obj.model.param.set('Study1D', x);

            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;
            figure(1)
            mphmesh(obj.model)
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.png')); % Change the path as needed
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.fig')); % Change the path as needed

            obj.model.study('std1').feature('time').set('tlist', 'range(0,1e-11,1.5e-9)');
            obj.model.study('std1').run;

            Wm = mphint2(obj.model, 'temw.Wm', 'line', 'selection', obj.Surf);
            Wmt = trapz(Wm);

            We = mphint2(obj.model, 'temw.We', 'line', 'selection', obj.Surf);
            Wet = trapz(We);

            Wm_in = mphint2(obj.model, 'temw.Wm', 'line', 'selection', obj.Surfin);
            Wm_int = trapz(Wm_in);

            We_in = mphint2(obj.model, 'temw.We', 'line', 'selection', obj.Surfin);
            We_int = trapz(We_in);        
        end


        function [Wt] = RunModel_Case(obj, Ip, x)
            mesh_min = obj.minmeshsize_nominal;
            mesh_max = obj.maxmeshsize_nominal;
            obj.model = obj.model;
            obj.Surf = obj.Surf;
            
            obj.model.param.set('maxmeshsize', mesh_max);
            obj.model.param.set('minmeshsize', mesh_min);
            obj.model.param.set('impact_location', Ip); % loaded as percentage over 1
            obj.model.param.set('Study1D', x);

            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;
            figure(1)
            if obj.plt
                mphmesh(obj.model)
            end
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.png')); % Change the path as needed
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.fig')); % Change the path as needed

            obj.model.study('std1').feature('time').set('tlist', 'range(0,1e-11,1.5e-9)');
            obj.model.study('std1').run;
            switch obj.objective
                case 'Wm'
                    Wm = mphint2(obj.model, 'temw.Wm', 'line', 'selection', obj.Surf);
                    Wt = trapz(Wm);
                case 'We'
                    We = mphint2(obj.model, 'temw.We', 'line', 'selection', obj.Surf);
                    Wt = trapz(We);
                case '(Wm+We)/Wm_in'
                    Wem = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surf);
                    W_in = mphint2(obj.model, 'temw.Wm', 'line', 'selection', obj.Surfin);
                    Wt = trapz(Wem./W_in);
                case '(Wm+We)/We_in'
                    Wem = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surf);
                    W_in = mphint2(obj.model, 'temw.We', 'line', 'selection', obj.Surfin);
                    Wt = trapz(Wem./W_in);
                case '(Wm+We)/(We_in+Wm_in)'
                    Wem = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surf);
                    W_in = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surfin);
                    Wt = trapz(Wem./W_in);
                case 'Wm+We'
                    Wem = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surf);
                    Wt = trapz(Wem);
                case '(Wm+We)/(We_in+Wm_in)*ymin'
                    Wem = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surf);
                    W_in = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surfin);
                    ymin = str2double(regexprep(char(obj.model.param.get('ymin')), '[^\d\.]', ''));
                    Wt = trapz(Wem./W_in) * (-ymin);
                otherwise
                    error('Unknown objective: %s', obj.objective);
            end

        end

        function saveData(obj, folder, mainName)
            % Save variables to .mat file
            matFileName = fullfile(folder, [mainName '_' obj.creationDate '.mat']);
            save(matFileName, 'obj');
            disp(['Variables saved to ', matFileName]);
            
            % Save variables to .csv file
            csvFileName = fullfile(folder, [mainName '_' obj.creationDate '.csv']);
            varNames = {'x', 'Object'};
            T = table(obj.x, obj.Object, 'VariableNames', varNames);
            writetable(T, csvFileName);
            disp(['Variables saved to ', csvFileName]);
        end


    end
end
