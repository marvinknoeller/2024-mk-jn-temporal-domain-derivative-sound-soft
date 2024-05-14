function reg = center_reg(ps,z,tvec)
%%%%%%%%%%%%%%
%% Input
%  ps: the spline parametrizing the boundary of the domain
%  z: the center point of the star domain
%  tvec: the values of t required to evaluate ps(t)
%% Output
%  reg: the regularization term
%%%%%%%%%%%%%%
psder = fnder(ps);
psval = fnval(ps,tvec/(2*pi)*ps.breaks(end));
psderval = fnval(psder, tvec/(2*pi)*ps.breaks(end));

A = 1/2*sum(2*pi/(length(tvec)) * (psval(1,:).*psderval(2,:) - psval(2,:).*psderval(1,:)));

int1 = 1/2*sum(2*pi/(length(tvec)) * ( psval(1,:).^2.*psderval(2,:)));
int2 = 1/2*sum(2*pi/(length(tvec)) * (-psval(2,:).^2.*psderval(1,:)));
center_ps = 1/A * [int1;int2];
reg = center_ps - z;

end