function penDer = center_regder(ps,hs,z,h,tvec)
%%%%%%%%%%%%%%
%% Input
%  ps: the spline parametrizing the boundary of the domain
%  hs: the spline parametrizing the pertubation of the boundary of the domain
%  z: the center point of the star domain (not needed!)
%  h: the (possible) shift of the center point
%  tvec: the values of t required to evaluate ps(t)
%% Output
%  penDer: the derivative of the regularization term
%%%%%%%%%%%%%%

psder = fnder(ps);
psval = fnval(ps,tvec/(2*pi)*ps.breaks(end));
psderval = fnval(psder, tvec/(2*pi)*psder.breaks(end));

hsder = fnder(hs);
hsval = fnval(hs,tvec/(2*pi)*hs.breaks(end));
hsderval = fnval(hsder,tvec/(2*pi)*hsder.breaks(end));


A = 1/2*sum(2*pi/(length(tvec)) * (psval(1,:).*psderval(2,:) - psval(2,:).*psderval(1,:)));
Ader = 1/2*sum( 2*pi/(length(tvec)) * (hsval(1,:).*psderval(2,:) - hsval(2,:).*psderval(1,:) ...
    + psval(1,:).*hsderval(2,:) - psval(2,:).*hsderval(1,:)) );


int1 = 1/2*sum(2*pi/(length(tvec)) * ( psval(1,:).^2.*psderval(2,:)));
int2 = 1/2*sum(2*pi/(length(tvec)) * (-psval(2,:).^2.*psderval(1,:)));

F1der = 1/2*sum( 2*pi/(length(tvec)) * (2*psval(1,:) .* hsval(1,:) .* psderval(2,:) ...
    + psval(1,:).^2 .* hsderval(2,:)) );
F2der = 1/2*sum( 2*pi/(length(tvec)) * (-2*psval(2,:) .* hsval(2,:) .* psderval(1,:) ...
    - psval(2,:).^2 .* hsderval(1,:)) );
% note that h = [0;0], when we derive w.r.t. radial components r.
% if we derive w.r.t. z, the h_1 = [1;0] and h_2 = [0;1] .
R1 = -1/A^2 * Ader * int1 + 1/A * F1der - h(1);
R2 = -1/A^2 * Ader * int2 + 1/A * F2der - h(2);
penDer = [R1;R2];

end