function [u,dnu]=pulseRK(g,x,x0,T,M,c)

%%%%%%%%%%%%%%
%% Input
%  g: function handle of 1d smooth function
%  x: points in space, where you want to have the wave
%  x0: sender / receiver
%  T: final time 
%  M: number of steps in time
%  c: coefficient from the RK method
%% Output
%  u: the pulse at the stages needed for the CQ at the points x
%%%%%%%%%%%%%%

tau = T/M;
m = length(c);
timepoints = get_time_points(T,tau, m, c);
u = zeros(size(x,1),m*M+1);

%%
diff = [];
for kk = 1:size(x0,2)
    t1 = x0(1,kk)-x(:,1);
    t2 = x0(2,kk)-x(:,2);
    diff = [diff;sqrt(t1.^2 + t2.^2)];
    ukk = zeros(size(x0,2)*size(x,1),m*M+1);
end
    parfor pp = 1 : length(diff) % for all points in space
        perm_ind1 = timepoints./diff(pp)>=1; % permitted indices
        evaluations1 = timepoints./diff(pp);
        evaluations1(~perm_ind1) = 0; % set non permitted to 0
        p1 = zeros(1,length(evaluations1));
        p1(perm_ind1) = arrayfun(@(t)...
        integral(@(theta)g(t*diff(pp)-diff(pp).*cosh(theta)),0,acosh(t))...
            ,evaluations1(perm_ind1));
        ukk(pp,:) = p1
    end
for kk = 1:size(x0,2)
    u = u + ukk((kk-1)*size(x,1)+1: kk*size(x,1),:);
end

end