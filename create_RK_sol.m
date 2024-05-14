function [uscat,lambdastages] = create_RK_sol(Var,Geometry)
%%%%%%%%%%%%%%
%% Input
%  Var, Geometry both cells containing the information about geometry,
%  parameters and incident wave
%% Output
%  uscat: the scattered wave at the points Geometry.observationpoints.
%         this might also be the domain derivative!
%  lambdastages: in the relevant situation, this will be up to some scaling
%                the Neumann data of the total wave at the boundary at all
%                the RK stages! This will be required for the boundary 
%                 condition of the domain derivative
%%%%%%%%%%%%%%
% Discrete geometry
g = starshapespline(Var.N,0,Geometry.ps);  % the usual grid
gp = starshapespline(Var.N,1/6,Geometry.ps);  % grid shifted by 1/6
gm = starshapespline(Var.N,-1/6,Geometry.ps); % grid shifted by -1/6
% construct the operators
V = CalderonCalculusHelmholtz(g,gp,gm); 
S = HelmholtzPotentials(g,Geometry.observationpoints);
%% Boundary integral equation coming from Single layer ansatz (Indirect)
if nargout == 2 % in this case we want to have lambdastages at the end
    d1   = Geometry.nobs;
    F = @(s,b)F_op(s,(-1)*S(s),V(s),b,d1,Geometry.mindist);
    [uscatstages ,lambdastages] = applyRKCQ(F,Geometry.beta0,d1,Var.T,Var.M,Var.A,Var.b,Var.c);
    m = length(Var.c);
    uscat = uscatstages(:,1:m:m*Var.M+1);
else
    d1   = Geometry.nobs;
    F = @(s,b)F_op(s,(-1)*S(s),V(s),b,d1,Geometry.mindist);
    uscatstages = applyRKCQ(F,Geometry.beta0,d1,Var.T,Var.M,Var.A,Var.b,Var.c);
    m = length(Var.c);
    uscat = uscatstages(:,1:m:m*Var.M+1,:);
    lambdastages = [];
end
end