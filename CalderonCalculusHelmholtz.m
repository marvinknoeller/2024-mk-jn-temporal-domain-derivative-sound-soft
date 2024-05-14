function [V] = CalderonCalculusHelmholtz(g,gp,gm)
%%%%%%%%%%%
%% Input:
%     g : geometry
%    gp : shifted geometry with epsilon = 1/6
%    gm : shifted geometry with epsilon = -1/6
%% Output:
%     V : averaged single layer operator V
%%%%%%%%%%%

Vp = CalderonCalculusHelmholtzHalf(g,gp);
Vm = CalderonCalculusHelmholtzHalf(g,gm);

Pp = 1/2; Pm = 1/2;

V =@(s) Pp*Vp(s) + Pm*Vm(s);

return

% Subfunction computing the two halves of the Calculus

function V = CalderonCalculusHelmholtzHalf(g,gp)
%%%%%%%%%%%
%% Input:
%     g : principal geometry
%    gp : companion geometry
%% Output:
%     V : single layer operator

DX = bsxfun(@minus,gp.midpt(:,1),g.midpt(:,1)');
DY = bsxfun(@minus,gp.midpt(:,2),g.midpt(:,2)');
D = sqrt(DX.^2+DY.^2);
V = @(s) 1i/4.*besselh(0,1,1i.*s.*D);


return


