clear all; close all; clc;
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% important parameters
video = 1;
Var = {};
[A,b,c] = RKdata(1);
Var.A = A;
Var.b = b;
Var.c = c;
reference_numpoints = 600;
reference_numtimes = 2000;
Var.T = 20;
Var.Tlag = 4;
Var.fork = 0;
val = 7;
Var.boxx = [-val -val val val];
Var.boxy = [val -val -val val];
% signal
HH = @mysignal;
HHp = @mysignalp;
amp1 = 3;
amp2 = 1;
signal = @(t) HH(amp1*t)- HH(amp2*t+2);
signalp= @(t) amp1*HHp(amp1*t) - amp2*HHp(amp2*t+2);
% save them
Var.signal = signal;
Var.signalp = signalp;
% direction of wave propogation
dir = 1/sqrt(2) * [1, 1];
Var.dir = dir;
% define the true scattering object
%% KITE
tt = linspace(0,2*pi-2*pi/1e5,1e5);
xx = cos(tt) + 0.65*cos(2*tt)-0.65;
yy = 1.5*sin(tt);
theta = pi/2;
res = [cos(theta) -sin(theta); sin(theta) cos(theta)]*[xx;yy];
xx = res(1,:); yy = res(2,:);
points = 1.5*[xx;yy];
points = [points, points(:,1)];
%%
Var.N = reference_numpoints;
Var.M = reference_numtimes;
ps = generatesplinehandles(points);
tline = linspace(0,2*pi-2*pi/Var.N,Var.N);
vals = fnval(ps,tline/(2*pi)*ps.breaks(end));
pexact = generatesplinehandles(points);
gp = starshapespline(Var.N,1/6,ps);
gm = starshapespline(Var.N,-1/6,ps);

% generate boundary data
Pp = 1/2; Pm = 1/2;
uincp = planewaveRK(signal,signalp,dir,Var.Tlag,gp.midpt,gp.normal,Var.T,Var.M,1,c);
uincm = planewaveRK(signal,signalp,dir,Var.Tlag,gm.midpt,gm.normal,Var.T,Var.M,1,c);
beta0 = (Pp*uincp+Pm*uincm);
% observation points
nobs = 8;
startobs = 0;
endobs =   2*pi;
tobs = linspace(startobs,endobs-endobs/nobs,nobs);
radius = 6;

X = radius * cos(tobs); X = X.';
Y = radius * sin(tobs); Y = Y.';

nobs = length(X);
Geometry.beta0 = beta0;
Geometry.ps = ps;
Geometry.observationpoints = [X,Y];
Geometry.nobs = nobs;
Geometry.mindist = minimaldistance(Var,Geometry);

uscat = direct_sound_soft_splineRK(Var,Geometry,video,'direct_vid'); fprintf("Direct problem: Done. \n")