function  PlotIterates(vals,ssvec,observationpoints,time_points,Name)
%%%%%%%%%%%%%%
%% Input
%  vals: the exact (unknown) boundary of the object
%  ssvec: tensor having the different iterates along the 3rd dimension
%  observationopints: the observationpoints
%  time_points: points in time
%  Name: String. How do you want to name the plots?
%% Output
%%%%%%%%%%%%%%
fig=figure(700);
fig.Position = [608 110 478 460];
% time_points = [1,3,6,11,21,size(ssvec,3)];
for ell = 1 : length(time_points)
        % figure
        % tline = linspace(0,2*pi,N);
        plot(vals(1,:), vals(2,:),'-b','LineWidth',2)
        hold on
        plot(ssvec(1,:,time_points(ell)), ssvec(2,:,time_points(ell)),'--r','LineWidth',2)
        ax = gca;
        ax.XLim = [-8,8];
        ax.YLim = [-8,8];
        ax.FontSize = 20;
        

        % s = plot(observationpoints(:,1),observationpoints(:,2),'dk','MarkerFaceColor','k','MarkerSize',5);
        plot(observationpoints(:,1),observationpoints(:,2)...
            ,'kd','MarkerFaceColor','w',...
            'MarkerSize',8,'LineWidth',2)
        if ell == 1
            legend("exact shape","initial guess",...
            'Interpreter','Latex','FontSize',24,...
            'Location','North','Orientation','horizontal')
        elseif ell == 6
            legend("exact shape",strcat("iterate ","$\partial D_{",num2str(time_points(ell)-1),"}$"),...
            'Interpreter','Latex','FontSize',24,...
            'Location','North','Orientation','horizontal')
        else
        legend("exact shape",strcat("iterate ","$\partial D_{",num2str(time_points(ell)-1),"}$"),...
            'Interpreter','Latex','FontSize',24,...
            'Location','North','Orientation','horizontal')
        end
        hold off
        ax.XTick = -7:2:7;
        ax.YTick = -7:2:7;
        drawnow
        saveas(gcf,strcat('plots/',Name,num2str(ell)),'epsc')
        % exportgraphics(ax,strcat('plots/',Name,num2str(ell),'.png'),'Resolution',300) 
end
end