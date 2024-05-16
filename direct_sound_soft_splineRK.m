function [uscat,lambdastages] = direct_sound_soft_splineRK(Var,Geometry,visualize,videoname)
%%%%%%%%%%%%%%
%% Input
%  Var, Geometry: Cells containing the information about geometry, waves
%  and parameters
%  visualize: Do you want to see a video? Fit the space and time parameters
%             as this might require a very long time else.
%  videoname: String defining the name of the video
%% Output
%  %  uscat: the scattered wave at the points Geometry.observationpoints.
%         this might also be the domain derivative!
%  lambdastages: in the relevant situation, this will be up to some scaling
%                the Neumann data of the total wave at the boundary at all
%                the RK stages! This will be required for the boundary 
%                 condition of the domain derivative
%%%%%%%%%%%%%%
if size(Geometry.beta0,1)==Var.N
    if nargout == 2
        [uscat,lambdastages] = create_RK_sol(Var,Geometry);
    else
        uscat = create_RK_sol(Var,Geometry);
        lambdastages = [];
    end
else
    uscat = create_RK_sol(Var,Geometry);
end


%% only go in here, if visualize == 1. Then, a video is generated that shows the time-dependent evolution
if nargin >=3
    if visualize >= 1
        Varvis = Var;
        Geometryvis = Geometry;
        [mesh,indoncurve] = create_triangulation_spline(Geometryvis.ps,Var);
        x = mesh.Nodes(1,:).';
        y = mesh.Nodes(2,:).';
        Xvis = x(:);
        Yvis = y(:);
        Geometryvis.observationpoints = [Xvis Yvis];
        Geometryvis.nobs = length(Xvis);
        Geometryvis.mindist = minimaldistance(Varvis,Geometryvis);
        uscatvis = create_RK_sol(Varvis,Geometryvis);
        % uincvis = planewaveRK(Var.signal,Var.signalp,Var.dir,Var.Tlag,[Xvis Yvis],0*[Xvis Yvis],Var.T,Var.M,0);
        if visualize == 1
            uincvis = planewaveRK(Var.signal,Var.signalp,Var.dir,Var.Tlag,[Xvis Yvis],0*[Xvis Yvis],Var.T,Var.M,1,Var.c);
        elseif visualize == 2
            x0 = Geometry.observationpoints.'; % Note that this is the same as x0!
            uincvis = pulseRK(Var.signal,[Xvis Yvis],x0,Var.T,Var.M,Var.c);
        else
            error('Choose either visualize = 1 or visualize = 2.');
        end
        m = length(Var.c);
        uincvis = uincvis(:,1:m:m*Var.M+1,:);
        utotalvis = zeros(length(Xvis),Var.M+1);
        utotalvis = uincvis+uscatvis;
        % This is a very simple post-processing step. We know that u=0 on the
        % boundary of the scattering object.
        [~,edges,~] = mesh.meshToPet;
        utotalvis(edges(1,indoncurve),:) = 0;
        % utotal = uinc+uscat;

        %% From here: Plot settings
        limu = max(utotalvis(:));
        % limd = min(utotalvis(:));
        Tri = mesh.Elements.';
        fig=figure(700);
        fig.Position = [608 110 478 460];
        if visualize == 1
            timepoints = [0, 4., 8.];
        elseif visualize == 2
            timepoints = [2.7, 5., 7.25];
        end
        for kk = 1 : length(timepoints)
            num = find(round(((1:Var.M+1)-1)*Var.T/Var.M,2)==timepoints(kk));

            trisurf(Tri(:,1:3),Xvis,Yvis,utotalvis(:,num));
            grid off
            ax = gca;
            ax.XLim = [-7,7];
            ax.YLim = [-7,7];
            ax.FontSize = 20;
            view(2)
            % ax.CLim = [-.5,.5];
            if visualize == 1
                ax.CLim = [-.7,.7];
            elseif visualize == 2
                ax.CLim = [-.3,.3];
            end
            colormap(bone)
            shading interp, shg
            if visualize == 1
            if timepoints(kk) == 0
                text(.9,5.4,limu+1,...
                    [strcat(" time t = ", num2str((num-1)*Var.T/Var.M)), '(Initial state)'],...
                    'FontSize',24,'Interpreter','Latex',...
                    'EdgeColor','w','LineStyle','-','LineWidth',1.5,'Margin',6,'Color','w')
            else
                text(.9,5.4,limu+1,...
                    strcat(" time t = ", num2str(round((num-1)*Var.T/Var.M,2),'%.2f')),...
                    'FontSize',24,'Interpreter','Latex',...
                    'EdgeColor','w','LineStyle','-','LineWidth',1.5,'Margin',6,'Color','w')
            end
            elseif visualize == 2
                if timepoints(kk) == 0
                text(.9,5.4,limu+1,...
                    [strcat(" time t = ", num2str((num-1)*Var.T/Var.M)), '(Initial state)'],...
                    'FontSize',24,'Interpreter','Latex',...
                    'EdgeColor','k','LineStyle','-','LineWidth',1.5,'Margin',6,'Color','k')
            else
                text(.9,5.4,limu+1,...
                    strcat(" time t = ", num2str(round((num-1)*Var.T/Var.M,2),'%.2f')),...
                    'FontSize',24,'Interpreter','Latex',...
                    'EdgeColor','k','LineStyle','-','LineWidth',1.5,'Margin',6,'Color','k')
            end
            end
            ax.XTick = -6:2:6;
            ax.YTick = -6:2:6;
            set(gcf,'renderer','painters') 
            saveas(gcf,strcat('plots/','PlotNr',num2str(kk),num2str(visualize)),'epsc')
            % exportgraphics(ax,strcat(videoname,'PlotNr',num2str(kk),'.png'),'Resolution',300)
            pause(1)
        end

    end
end

end