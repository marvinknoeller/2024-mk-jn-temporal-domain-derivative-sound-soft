function [mesh,indoncurve] = create_triangulation_spline(ps,Var)
%% this is just for the visualization, i.e. for generating a video!
%%%%%%%%%%%%%%
%% Input
%  ps: the spline parametrizing the boundary of the domain
%  Var: the cell containing the variables
%% Output
%  mesh: a mesh around the scattering object
%  indoncurve: indices that define the points on the boundary of the object
%              This will be needed for a simple postprocessing step
%%%%%%%%%%%%%%
N = Var.N;
boxx = Var.boxx;
boxy = Var.boxy;
tline = linspace(0,2*pi-2*pi/N,N);
ss = fnval(ps,tline/(2*pi)*ps.breaks(end));
psder = fnder(ps);
ssder = fnval(psder,tline/(2*pi)*ps.breaks(end));
normals = [ssder(2,:); -ssder(1,:)];
unitnormals = normals./sqrt(sum(normals.^2,1));
pgon = polyshape({boxx, ss(1,:)+1e-2*unitnormals(1,:)}, ...
    {boxy, ss(2,:)+1e-2*unitnormals(2,:)});

tr = triangulation(pgon);
model = createpde;
tnodes = tr.Points';
telements = tr.ConnectivityList';

geometryFromMesh(model,tnodes,telements);
tolob = 0.05;
mesh = generateMesh(model,'Hmax',.07,...
    "Hedge",{5,tolob},...
    'GeometricOrder','linear');
[~,edges,~] = mesh.meshToPet;

maximal_curve_segments = max(edges(5,:));
indoncurve = [];
for kk = 5:maximal_curve_segments
    indoncurve = [indoncurve, find(edges(5,:) == kk)]; %find all the indices on the curve
end

end