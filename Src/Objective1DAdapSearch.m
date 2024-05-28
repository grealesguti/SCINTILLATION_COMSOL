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
        savename
        counter
        SimpStol
        I0
        Iend
    end
    
    methods
        function obj = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surfin, model, plt, objective, savename,SimpStol,I0,Iend)
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
            obj.creationDate = datestr(now, 'yyyymmdd_HHMM');
            obj.objective = objective;
            obj.savename = savename;
            obj.counter=1;
            obj.SimpStol=SimpStol;
            obj.I0=I0;
            obj.Iend=Iend;
        end
        
        function WI = compute(obj, x)
            % Print vector x to the screen
            fprintf('Vector x: ');
            disp(x);

            LYSO_L = 0.0275;

            x0 = 0.05;
            xend = 1.95;

            if (x0 < x) && (x < xend)
                FUN = @(Ip) obj.ObjectiveQuad_1D(Ip, x);
                [WI, ~, ~, t, y] = adaptiveSimpson(FUN, obj.I0, obj.Iend, 'parts', 2, 'tol', obj.SimpStol);
                name_ty=append('Rst/', obj.savename, '_x_',num2str(x, '%.3f'),'_date_', obj.creationDate);
                saveData(name_ty, 't', t, 'y', y);
                disp(['t & y Variables saved to ', name_ty]);
                WI = -WI;
                disp(['Objective value: ', sprintf('%.2e', WI)]);
            else
                WI = 1000;
            end
        end


        function WI = compute_Nvar(obj, x, Nvar)
            % Print vector x to the screen
            fprintf('Vector x: ');
            disp(x);

            LYSO_L = 0.0275;
            x0 = 0.05;
            xend = 1.95;

                FUN = @(Ip) obj.ObjectiveQuad_Nvar(Ip, x, Nvar);
                [WI, ~, ~, t, y] = adaptiveSimpson(FUN, obj.I0, obj.Iend, 'parts', 2, 'tol', obj.SimpStol);
                if length(x)==1
                    name_ty=append('Rst/', obj.savename, '_x_',num2str(x, '%.3f'),'_date_', obj.creationDate);
                    saveData(name_ty, 't', t, 'y', y);
                else 
                    count=obj.counter;
                    obj.counter=obj.counter+1;
                    name_ty=append('Rst/', obj.savename, '_x_',num2str(obj.counter),'_date_', obj.creationDate);
                    saveData(name_ty, 't', t, 'y', y);
                end
                disp(['t & y Variables saved to ', name_ty]);
                WI = -WI;
                disp(['Objective value: ', sprintf('%.2e', WI)]);

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

        function WI = ObjectiveQuad_Nvar(obj, Ip, x, Nvar)
            % Initialize Run function
            FUN = @(xx) obj.RunModel_Case_Nvar(xx, x, Nvar);
            WI = zeros(length(Ip), 1);
            for i = 1:length(Ip)
                disp(['Solving impact number: ', num2str(i), ' Out of: ', num2str(length(Ip)), '. Location: ', num2str(Ip(i))]);
                WI(i) = FUN(Ip(i));
            end
            WI = WI';
        end        


        function [Wt,other] = RunModel_Case_Nvar(obj, Ip, x, Nvar)
            mesh_min = obj.minmeshsize_nominal;
            mesh_max = obj.maxmeshsize_nominal;
            obj.model = obj.model;
            obj.Surf = obj.Surf;
            other=[];
            
            obj.model.param.set('maxmeshsize', mesh_max);
            obj.model.param.set('minmeshsize', mesh_min);
            obj.model.param.set('impact_location', Ip); % loaded as percentage over 1

            for i=1:Nvar/2
                obj.model.param.set(append('p',num2str(i-1)), x(i));
                obj.model.param.set(append('m',num2str(i-1)), x(i+Nvar/2));
            end

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
                    other=We;
                case 'We/We_in'
                    Wem = mphint2(obj.model, 'temw.We', 'line', 'selection', obj.Surf);
                    W_in = mphint2(obj.model, 'temw.We', 'line', 'selection', obj.Surfin);
                    Wt = trapz(Wem./W_in);
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

        function [Wt,other] = RunModel_Case(obj, Ip, x)
            mesh_min = obj.minmeshsize_nominal;
            mesh_max = obj.maxmeshsize_nominal;
            obj.model = obj.model;
            obj.Surf = obj.Surf;
            other=[];

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
                    other=We;
                case '(Wm+We)/Wm_in'
                    Wem = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surf);
                    W_in = mphint2(obj.model, 'temw.Wm', 'line', 'selection', obj.Surfin);
                    Wt = trapz(Wem./W_in);
                case 'Win'
                    W_in = mphint2(obj.model, 'temw.Wm+temw.We', 'line', 'selection', obj.Surfin);
                    Wt = trapz(W_in);
                case 'We_in'
                    W_in = mphint2(obj.model, 'temw.We', 'line', 'selection', obj.Surfin);
                    Wt = trapz(W_in);
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

        function [c, ceq] = ComsolVolumeConstraint_ga(obj, x, MaxVal, Nvar)
            for i = 1:Nvar/2
                obj.model.param.set(append('p', num2str(i-1)), x(i));
                obj.model.param.set(append('m', num2str(i-1)), x(i+Nvar/2));
            end
        
            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;
        
            obj.model.component('comp1').geom('geom1').measure().selection().init(2);
            obj.model.component('comp1').geom('geom1').measure().selection().set('uni1',[1,2]);
            VarComsolVal = obj.model.component('comp1').geom('geom1').measure().getVolume();
        
            % Compute constraint values
            c = VarComsolVal/MaxVal - 1; % Inequality constraint
            ceq = []; % No equality constraints
        end


        function [ConVal]=ComsolVolumeConstraint(obj, x, MaxVal, Nvar)

            for i=1:Nvar/2
                obj.model.param.set(append('p',num2str(i-1)), x(i));
                obj.model.param.set(append('m',num2str(i-1)), x(i+Nvar/2));
            end

            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;

            obj.model.component('comp1').geom('geom1').measure().selection().init(2);
            obj.model.component('comp1').geom('geom1').measure().selection().set('uni1',[1,2]);
            VarComsolVal=obj.model.component('comp1').geom('geom1').measure().getVolume();
            %VarComsolVal=obj.model.component('comp1').geom(Domains).measure().getVolume();
            %VarComsolVal=obj.model.param.get(ComsolVarName); % I KNOW THIS IS WRONG
            ConVal = VarComsolVal/MaxVal-1;

            fprintf('Vol Value: %d, Constraint: %d\n',VarComsolVal, ConVal)
        end

        function [Distance]=ComsolDistancePts(obj,pto1,pto2)
            % RUN
            obj.model.component('comp1').geom('geom1').run;
            geom = obj.model.geom('geom1');
            vertex = geom.getVertexCoord;
            pto1_coord = vertex(:,pto1);
            pto2_coord = vertex(:,pto2);
            Distance = abs(sqrt(sum((pto1_coord - pto2_coord).^2)));
            %obj.model.component('comp1').geom('geom1').measure().selection().init(2);
            %obj.model.component('comp1').geom('geom1').measure().selection().set('dif2',[1]);
            %VarComsolVal=obj.model.component('comp1').geom('geom1').measure().getVolume();
            %VarComsolVal=obj.model.component('comp1').geom(Domains).measure().getVolume();
            %VarComsolVal=obj.model.param.get(ComsolVarName); % I KNOW THIS IS WRONG
            %ConVal = VarComsolVal/MaxVal-1;

            fprintf('Distance Value: %d\n',Distance)
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
