function  PlotIteratesTotal(vals,ssvec,observationpoints,time_point,Name,leg)
%%%%%%%%%%%%%%
%% Input
%  vals: the exact (unknown) boundary of the object
%  ssvec: tensor having the different iterates along the 3rd dimension
%  observationopints: the observationpoints
%  time_points: points in time
%  Name: String. How do you want to name the plots?
%  leg: Do you want a legend?
%% Output
%%%%%%%%%%%%%%
fig=figure(700);
fig.Position = [608 110 478 460];
% time_points = [1,3,6,11,21,size(ssvec,3)];
        % figure
        % tline = linspace(0,2*pi,N);
        s1 = plot(ssvec(1,:,1), ssvec(2,:,1),'-.k','LineWidth',2);
        hold on
        s2 = plot(vals(1,:), vals(2,:),'-b','LineWidth',2);
        s3 = plot(ssvec(1,:,time_point), ssvec(2,:,time_point),'--r','LineWidth',2);
        ax = gca;
        ax.XLim = [-7,7];
        ax.YLim = [-7,7];
        ax.FontSize = 20;
        

        s4 = plot(observationpoints(:,1),observationpoints(:,2),'ok','MarkerFaceColor','k','MarkerSize',5);
        
        % if ell == 1
        %     legend("exact shape","initial guess",...
        %     "$\mathbf{z}_j,\; j=1,\dots,8$",'Interpreter','Latex','FontSize',20,...
        %     'Location','NorthWest')
        % elseif ell == 6
            % if  leg == 1
            %     legend([s2,s3,s4,s1],"exact shape",strcat("final approx. ","$\partial D_{",num2str(time_point-1),"}$"),...
            %     "$\mathbf{z}_j,\; j=1,\dots,8$",'initial guess','Interpreter','Latex','FontSize',20,...
            %     'Location','NorthWest')
            % end
            if  leg == 1
                legend([s2,s3,s1],"exact shape",strcat("final approx. ","$\partial D_{",num2str(time_point-1),"}$"),...
                'initial guess','Interpreter','Latex','FontSize',25,...
                'Location','NorthWest')
            end
        % else
        % legend("exact shape",strcat("iterate $\partial D_{",num2str(time_point-1),"}$"),...
        %     "$\mathbf{z}_j,\; j=1,\dots,8$",'Interpreter','Latex','FontSize',20,...
        %     'Location','NorthWest')
        % end
        
        hold off
        ax.XTick = -6:2:6;
        ax.YTick = -6:2:6;
        drawnow
        saveas(gcf,strcat('plots/',Name),'epsc')
        % exportgraphics(ax,strcat('plots/',Name,'.png'),'Resolution',300) 
end