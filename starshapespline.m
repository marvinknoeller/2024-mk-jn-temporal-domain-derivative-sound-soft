function  g = starshapespline(N,ep,ps)
%%%%%%%%%%%%%%
%% Input
%   N: discretization parameter
%   ep: epsilon parameter
%   ps: the spline parametrizing the boundary of the object
%% Output
%   g: the discrete geometry
%%%%%%%%%%%%%%

h = 1/N;
t = h*(0:N-1); t = t+ep*h; t=t(:);
tau  = 2*pi*(t-0.5*h);
t    = 2*pi*t;
t = mod(t,2*pi);
tau = mod(tau,2*pi);
g.midpt = fnval(ps,t/(2*pi)*ps.breaks(end)).';
g.brkpt = fnval(ps,tau/(2*pi)*ps.breaks(end)).';

dpp = fnder(ps);

xp = fnval(dpp,t/(2*pi)*ps.breaks(end)).';
g.normal  = 2*pi*h*[xp(:,2) -xp(:,1)];
g.unitnormal  = g.normal./sqrt(sum(abs(g.normal).^2,2));
g.next    = [2:N 1];
g.comp    = [1];

return