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
        intshape
        deltaY
        ImpactObject
        Weimag
        ydata
        ymin
        We
        mindist
        parts
    end
    
    methods
        function obj = Objective1DAdapSearch(minmeshsize_nominal, maxmeshsize_nominal, Surf, Surfin, model, plt, objective, savename,SimpStol,I0,Iend,int,deltaY,IObj, Weimag)
            import com.comsol.model.util.* 
            obj.plt = plt;
            if plt
                ModelUtil.showProgress(true);
            end
            obj.parts = 2;
            obj.minmeshsize_nominal = minmeshsize_nominal;
            obj.maxmeshsize_nominal = maxmeshsize_nominal;
            obj.Surf = Surf;
            obj.Surfin = Surfin;        
            fprintf('Opening MPH')
            obj.model = mphopen(model);
            fprintf('Opened')
            obj.creationDate = datestr(now, 'yyyymmdd_HHMM');
            obj.objective = objective;
            obj.savename = savename;
            obj.counter=1;
            obj.SimpStol=SimpStol;
            obj.I0=I0;
            obj.Iend=Iend;
            if strcmp(int,'2D')
                obj.intshape = 'surface';
            else
                obj.intshape = 'line';
            end
            obj.deltaY=deltaY;
            obj.ImpactObject=IObj;
            obj.Weimag=Weimag;
            obj.ymin =4.352009e-02/1.5;
            obj.ydata =[];
            obj.mindist = 2/3;
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
                [WI, ~, ~, t, y] = adaptiveSimpson(FUN, obj.I0, obj.Iend, 'parts', obj.parts, 'tol', obj.SimpStol);
                name_ty=append('Rst/', obj.savename, '_x_',num2str(x, '%.3f'),'_date_', obj.creationDate);
                if obj.plt
%                     figure(2)
%                     hold on 
%                     plot(t,y)
%                     figure(1)
                    % Sort t in increasing order and reorder y accordingly
                    [t_sorted, idx] = sort(t);
                    y_sorted = y(idx);
                
                    figure(2)
                    hold on 
                    plot(t_sorted, y_sorted)
                
                    figure(1)                    
                end
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
            if isempty(obj.We)
                FUN = @(Ip) obj.ObjectiveQuad_Nvar(Ip, x, Nvar);
                [WI, ~, ~, t, y] = adaptiveSimpson(FUN, obj.I0, obj.Iend, 'parts', obj.parts, 'tol', obj.SimpStol);
                if length(x)==1
                    name_ty=append('Rst/', obj.savename, '_x_',num2str(x, '%.3f'),'_date_', obj.creationDate);
                    saveData(name_ty, 't', t, 'y', y);
                    obj.ydata = y;
                else 
                    count=obj.counter;
                    obj.counter=obj.counter+1;
                    name_ty=append('Rst/', obj.savename, '_x_',num2str(obj.counter),'_date_', obj.creationDate);
                    saveData(name_ty, 't', t, 'y', y);
                    obj.ydata = y;
                end
                disp(['t & y Variables saved to ', name_ty]);
                WI = -WI;
                disp(['Objective value: ', sprintf('%.2e', WI)]);
            else
                WI=obj.We;
                disp(['Objective value: ', sprintf('%.2e', WI)]);
                obj.We=[];
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

        function R2 = ObjectiveRhoCond(obj, Ip, x)
            % Initialize Run function
            FUN = @(xx) obj.RunModelRhoCond_Case(xx, x);
            WI = zeros(length(Ip), 1);
            for i = 1:length(Ip)
                disp(['Solving impact number: ', num2str(i), ' Out of: ', num2str(length(Ip)), '. Location: ', num2str(Ip(i))]);
                WI(i) = FUN(Ip(i));
            end
            G4rst= obj.funG4RhoCond(Ip);
            error =  WI./WI(1) - G4rst';
            R2 = sum(error.^2)*100;
            figure(2)
            hold on
            plot(WI./WI(1))
            plot(G4rst)
            hold off
            disp(['Error R2: ', num2str(R2)]);
            disp(error)

        end        


        function R2 = ObjectiveRhoCondIncr(obj,fig, Ip, Iv, rho)
            % Initialize Run function
            FUN = @(xx) obj.RunModelgen(xx);
            obj.model.param.set('Study1D',Iv)
            obj.model.param.set('cond_LYSO',rho)
            WI = zeros(length(Ip), 1);

            for i = 1:length(Ip)
                disp(['Solving impact number: ', num2str(i), ' Out of: ', num2str(length(Ip)), '. Location: ', num2str(Ip(i))]);
                WI(i) = FUN(Ip(i));
            end
            figure(fig)
            hold on
            plot(WI)
            %plot(G4rst)
            hold off
        end        



        function [G4norm] = funG4RhoCond(obj,Ip)
            G4norm = 1.0003-0.0075556*Ip+0.065409*Ip.^2;
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

            %obj.model.study('std1').feature('time').set('tlist', 'range(0,1e-11,1.5e-9)');

            % Start the timer
            tic;
            obj.model.study('std1').run;
            elapsedTime = toc;
            fprintf('Time taken to run: %.2f seconds\n', elapsedTime);
            
            We = mphint2(obj.model, obj.objective, obj.intshape, 'selection', obj.Surf);
                    if length(We)>1
                        Wt = trapz(We);
                    else
                        Wt=We;
                    end
                    if not(isempty(obj.Weimag))
                        Weimagval = mphint2(obj.model, 'imag(temw.We)', 'surface', 'selection', obj.Weimag);
                        Wt=Wt-abs(trapz(Weimagval));
                    end
                    if obj.deltaY>0
                        obj.model.component('comp1').geom('geom1').measure().selection().init(1);
                        obj.model.component('comp1').geom('geom1').measure().selection().set(obj.ImpactObject,obj.deltaY);
                        Yval=obj.model.component('comp1').geom('geom1').measure().getVolume();
                        Wt=Wt/(Yval);
                    end

                    fprintf('The value of Wt is:  %e\n', Wt);
                    other=We;

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

            tic;
            obj.model.study('std1').run;
            elapsedTime = toc;
            fprintf('Time taken to run: %.2f seconds\n', elapsedTime);

                    We = mphint2(obj.model, obj.objective, obj.intshape, 'selection', obj.Surf);
                    if length(We)>1
                        Wt = trapz(We);
                    else
                        Wt=We;
                    end
                    if not(isempty(obj.Weimag))
                        Weimagval = mphint2(obj.model, 'imag(temw.We)', 'surface', 'selection', obj.Weimag);
                        Weimagvaltrap=trapz(Weimagval);
                        Wt=Wt-abs(Weimagvaltrap);
                    end
                    if obj.deltaY==-1
                        Wt=Wt/(x+2*(1-x)*Ip);
                    elseif obj.deltaY>0
                        obj.model.component('comp1').geom('geom1').measure().selection().init(1);
                        obj.model.component('comp1').geom('geom1').measure().selection().set(obj.ImpactObject,obj.deltaY);
                        Yval=obj.model.component('comp1').geom('geom1').measure().getVolume();
                        Wt=Wt/(Yval);
                    end
                    fprintf('The value of Wt is:  %e\n', Wt);
                    other=We;

        end
%%
        function [Wt,other] = RunModelRhoCond_Case(obj, Ip, x)
            mesh_min = obj.minmeshsize_nominal;
            mesh_max = obj.maxmeshsize_nominal;
            obj.model = obj.model;
            obj.Surf = obj.Surf;
            other=[];

            obj.model.param.set('maxmeshsize', mesh_max);
            obj.model.param.set('minmeshsize', mesh_min);
            obj.model.param.set('impact_location', Ip); % loaded as percentage over 1
            obj.model.param.set('Study1D', 1);
            obj.model.param.set('cond_LYSO', x);

            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;
            figure(1)
            if obj.plt
                mphmesh(obj.model)
            end
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.png')); % Change the path as needed
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.fig')); % Change the path as needed

            %obj.model.study('std1').feature('time').set('tlist', 'range(0,1e-11,1.5e-9)');
            tic;
            obj.model.study('std1').run;
            elapsedTime = toc;
            fprintf('Time taken to run: %.2f seconds\n', elapsedTime);
            
            
            We = mphint2(obj.model, obj.objective, obj.intshape, 'selection', obj.Surf);
                    Wt = trapz(We);
                    fprintf('The value of Wt is:  %e\n', Wt);
                    other=We;

        end



        function [Wt,other] = RunModelgen(obj, Ip)
            mesh_min = obj.minmeshsize_nominal;
            mesh_max = obj.maxmeshsize_nominal;
            obj.model = obj.model;
            obj.Surf = obj.Surf;
            other=[];

            obj.model.param.set('maxmeshsize', mesh_max);
            obj.model.param.set('minmeshsize', mesh_min);
            obj.model.param.set('impact_location', Ip); % loaded as percentage over 1

            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;
            figure(1)
            if obj.plt
                mphmesh(obj.model)
            end
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.png')); % Change the path as needed
            %saveas(gcf, append('./MeshConv/mesh_',num2str(Ip),'_',num2str(mesh_max),'_',num2str(mesh_min),'.fig')); % Change the path as needed

                        tic;
            obj.model.study('std1').run;
            elapsedTime = toc;
            fprintf('Time taken to run: %.2f seconds\n', elapsedTime);

                      We = mphint2(obj.model, obj.objective, obj.intshape, 'selection', obj.Surf);
                    Wt = trapz(We);
                    fprintf('The value of Wt is:  %e\n', Wt);
                    other=We;

        end

        %%

        function [c, ceq] = ComsolVolumeConstraint_ga(obj, x, MaxVal, Nvar)
            for i = 1:Nvar/2
                obj.model.param.set(append('p', num2str(i-1)), x(i));
                obj.model.param.set(append('m', num2str(i-1)), x(i+Nvar/2));
            end
            ceq = [];
            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;
        
            obj.model.component('comp1').geom('geom1').measure().selection().init(2);
            obj.model.component('comp1').geom('geom1').measure().selection().set('uni1',[1,2]);
            VarComsolVal = obj.model.component('comp1').geom('geom1').measure().getVolume();
        
            % Compute constraint values
            cv = VarComsolVal/MaxVal - 1; % Inequality constraint
            if isempty(obj.ymin)
                c=zeros(length(x)/2+1);
            else
                c=zeros(length(x)/2+2);
                c(end)= obj.ydata/obj.ymin-1;
            end
           
            c(1)=cv;
            for i = 1:length(x)/2
               c(i+1) = (x(i)+x(i+length(x)/2))/obj.mindist-1;
            end
            ceq = []; % No equality constraints
        end


        function [c, ceq]=ComsolVolumeConstraint(obj, x, MaxVal, Nvar)
            ceq = [];
            for i=1:Nvar/2
                obj.model.param.set(append('p',num2str(i-1)), x(i));
                obj.model.param.set(append('m',num2str(i-1)), x(i+Nvar/2));
            end

            if isempty(obj.We)
                FUN = @(Ip) obj.ObjectiveQuad_Nvar(Ip, x, Nvar);
                [WI, ~, ~, t, y] = adaptiveSimpson(FUN, obj.I0, obj.Iend, 'parts', obj.parts, 'tol', obj.SimpStol);
                obj.We=WI;
            end

            % RUN
            obj.model.component('comp1').geom('geom1').run;
            obj.model.component('comp1').mesh('mesh1').run;

            obj.model.component('comp1').geom('geom1').measure().selection().init(2);
            obj.model.component('comp1').geom('geom1').measure().selection().set('dif2',[1,2]);
            VarComsolVal=obj.model.component('comp1').geom('geom1').measure().getVolume();
            %VarComsolVal=obj.model.component('comp1').geom(Domains).measure().getVolume();
            %VarComsolVal=obj.model.param.get(ComsolVarName); % I KNOW THIS IS WRONG
            cv = VarComsolVal/MaxVal-1;


            % Compute constraint values
            if isempty(obj.ymin)
                c=zeros(length(x)/2+1,1);
            else
                c=zeros(length(x)/2+2,1);
                c(end)= -min(y)/obj.ymin+1;
            end
           
            c(1)=cv;
            for i = 1:length(x)/2
               c(i+1) = -(x(i)+x(i+length(x)/2))/obj.mindist+1;
            end

            fprintf('Vol Value: %d, Constraint: %d\n',VarComsolVal, cv)
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
