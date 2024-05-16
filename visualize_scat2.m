clear all; close all; clc;
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% important parameters
video = 2;
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
signal = @(t) HH(2*(t-1));
% save them
Var.signal = signal;
% define the true scattering object
load('pidgy.mat');
points = 2.5*[points, points(:,1)];
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
x0 = [-2.5, 2.5, -2.5, 2.5, 5;...
    -5.5, -5.5, 5.5, 5.5, 0];
tic
uincp = pulseRK(signal,gp.midpt,x0,Var.T,Var.M,c);
uincm = pulseRK(signal,gm.midpt,x0,Var.T,Var.M,c);
toc
beta0 = (Pp*uincp+Pm*uincm);
% observation points
X = x0(1,:).';
Y = x0(2,:).';
nobs = length(X);
Geometry.beta0 = beta0;
Geometry.ps = ps;
Geometry.observationpoints = [X,Y];
Geometry.nobs = nobs;
Geometry.mindist = minimaldistance(Var,Geometry);

uscat = direct_sound_soft_splineRK(Var,Geometry,video,'direct_vid'); fprintf("Direct problem: Done. \n")