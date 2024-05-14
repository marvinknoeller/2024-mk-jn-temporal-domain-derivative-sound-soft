function pen = curvepenalty(ps,tvec)
%%%%%%%%%%%%%%
%% Input
%  ps: the spline parametrizing the boundary of the domain
%  tvec: the values of t required to evaluate ps(t)
%% Output
%  pen: the regularization term
%%%%%%%%%%%%%%
psder = fnder(ps);
psderder = fnder(psder);

psderval = fnval(psder, tvec);
psderderval = fnval(psderder, tvec);

kappa1 = psderval(1,:).*psderderval(2,:) - psderderval(1,:).*psderval(2,:);

normpder = sqrt(sum(psderval.^2,1));

pen = sqrt(2*pi/(length(tvec))) * kappa1./ normpder.^(5/2);


end