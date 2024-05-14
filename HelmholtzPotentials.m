function [SL]=HelmholtzPotentials(g,z)
%%%%%%%%%%%%%%
%% Input
%  g: the geometry of the boundary of the scatterer
%  z: these are the output points
%% Output
%  SL: the function handle of the single layer potential
%%%%%%%%%%%%%%
RX = bsxfun(@minus,z(:,1),g.midpt(:,1)');   
RY = bsxfun(@minus,z(:,2),g.midpt(:,2)');
R = sqrt(RX.^2+RY.^2);  
SL = @(s) 1i/4*besselh(0,1,1i*s*R);

return