function mindist = minimaldistance(Var,Geometry)
%%%%%%%%%%%%%%
%% Input
%  Var, Geometry : Cells containing the information about geometry, waves
%  and parameters
%% Output
%  mindist : the minimal distance from the receivers to the obstacle
%%%%%%%%%%%%%%
g = starshapespline(Var.N,0,Geometry.ps);
points = g.midpt;
obs = Geometry.observationpoints;
for kk = 1 : length(obs)
    dist(kk) = min(sqrt(sum(abs(points - obs(kk,:)).^2,2)));
end
mindist = min(dist);
end