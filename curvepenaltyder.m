function derpen = curvepenaltyder(ps,hs,tvec)
%%%%%%%%%%%%%%
%% Input
%  ps: the spline parametrizing the boundary of the domain
%  hs: the spline parametrizing the perturbation of the boundary of the domain
%  tvec: the values of t required to evaluate ps(t)
%% Output
%  derpen: derivative of the penalty term
%%%%%%%%%%%%%%
psder = fnder(ps);
psderder = fnder(psder);

psderval = fnval(psder, tvec);
psderderval = fnval(psderder, tvec);

kappa1 = psderval(1,:).*psderderval(2,:) - psderderval(1,:).*psderval(2,:);
normpder = sqrt(sum(psderval.^2,1));


hsder = fnder(hs);
hsderder = fnder(hsder);

hsderval = fnval(hsder, tvec);
hsderderval = fnval(hsderder, tvec);

kappa1der = psderval(1,:).*hsderderval(2,:) + hsderval(1,:).*psderderval(2,:) ...
    - psderderval(1,:).*hsderval(2,:) - hsderderval(1,:).*psderval(2,:);

derpen = sqrt(2*pi/(length(tvec))) * ( kappa1der./(normpder.^(5/2)) - 5/2*kappa1./(normpder.^(9/2)).*sum(psderval.*hsderval,1) );

end