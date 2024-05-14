clear all; close all; clc;
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
parpool(64)
for it = 1 : 5
    rads = .5;
    save_results = 1;
    name = strcat("Paper11_",num2str(it));
    %% important parameters
    video = 0;
    regpar = 2e-2;
    cenpar = 5e-1; 
    degcoeff = 40;
    
    golden_section = 0;
    if golden_section == 1
        numsplits = 5;
    end
    if it >= 2
        noiselevel = .3; % 20%
    end


    %% direct problem first.
    Var = {};
    [A,b,c] = RKdata(1);
    Var.A = A;
    Var.b = b;
    Var.c = c;
    reference_numpoints = 2000;
    reference_numtimes = 20000;
    numpoints = 1000;
    numtimes = 5000;
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
    uscat = uscat(:,1:reference_numtimes/numtimes:end);
    Var.N = numpoints;
    Var.M = numtimes;
    
    Var.tau = Var.T/Var.M;
    % Add additional noise to the reference solution
    if it >=2
        rng(123) % set the seed.
        normuscat = errorongivenpoints(uscat,Var.tau);
        noise = -1 + 2*rand(size(uscat));
        noise = noise/errorongivenpoints(noise,Var.tau);
        uscat = uscat + noise*normuscat*noiselevel;
    end
    %% inverse problem
    % 1. step: Generate the initial grid
    % parametrize the radius
    tline = linspace(0,2*pi-2*pi/Var.N,Var.N);
    part = linspace(0,2*pi-2*pi/degcoeff,degcoeff);
    rvecell = rads*ones(1,degcoeff);
    % z_0 = [0;0.];
    z_0 = [-1;-1.5];
    pointsell = rvecell.*[cos(part);sin(part)] + z_0;
    pointsell = [pointsell, pointsell(:,1)];
    psell = generatesplinehandles(pointsell);
    ss = fnval(psell,tline/(2*pi)*psell.breaks(end));
    ssvec(:,:,1) = ss;
    g = starshapespline(Var.N,0,psell);
    gp = starshapespline(Var.N,1/6,psell);
    gm = starshapespline(Var.N,-1/6,psell);
    % update the Geometry
    Geometry.ps = psell;
    Geometry.mindist = minimaldistance(Var,Geometry);
    % 2. step : Solve the forward problem
    % generate boundary data
    Pp = 1/2; Pm = 1/2;
    uincp = planewaveRK(signal,signalp,dir,Var.Tlag,gp.midpt,gp.normal,Var.T,Var.M,1,c);
    uincm = planewaveRK(signal,signalp,dir,Var.Tlag,gm.midpt,gm.normal,Var.T,Var.M,1,c);
    beta0 = (Pp*uincp+Pm*uincm);
    % update the function on the boundary
    Geometry.beta0 = beta0;
    [uscat_ell,lambdastages] = direct_sound_soft_splineRK(Var,Geometry);
    uscatVar{1} = uscat_ell;
    
    psder = fnder(psell);
    psderval = fnval(psder, tline/(2*pi)*psder.breaks(end));
    unu_ell =  Var.N/(2*pi)*(1./sqrt(sum(psderval.^2,1))).'.*lambdastages;
    
    
    % 3. step: updating. Perform the Gauss-Newton step

    %%
    j = 1;
    relative_movement = 1;
    tol = 1e-2;
    %%
    while j < 5000
        if relative_movement < tol
           break
        end
        j
        error_ell = errorongivenpoints(uscat_ell - uscat,Var.tau) / errorongivenpoints(uscat,Var.tau); %this is not the functional that is supposed to be minimized - but the sqrt
        fprintf("perform a Newton step\n")
        fprintf(strcat("Error = ", num2str(error_ell), "\n"))
        boundarymat = [];
        for jj = 1 : degcoeff   % move radius and keep the center fixed
            rvecder = zeros(1,degcoeff);
            rvecder(jj) = 1;
            pointsder = rvecder.*[cos(part);sin(part)];
            pointsder = [pointsder, pointsder(:,1)];
            pselljj = generatesplinehandles(pointsder);

            ss = fnval(pselljj,tline/(2*pi)*pselljj.breaks(end));
            hnu = sum( ss .* g.unitnormal.',1);
            boundarydata = hnu.' .* unu_ell;
            %%
            boundarymat = [boundarymat; boundarydata];
            regprime(:,jj) = curvepenaltyder(psell,pselljj,tline);
            centerprime(:,jj) = center_regder(psell,pselljj,z_0,[0;0],tline);
        end

        for jj2 = 1 : 2 %x and y component
            if jj2 == 1
                pointsder = [ones(1,length(part)) ; zeros(1,length(part))];
            else
                pointsder = [ zeros(1,length(part)); ones(1,length(part))];
            end
            pointsder = [pointsder, pointsder(:,1)];
            pselljj = generatesplinehandles(pointsder);
            ss = fnval(pselljj,tline/(2*pi)*pselljj.breaks(end));
            hnu = sum( ss .* g.unitnormal.',1);
            boundarydata = hnu.' .* unu_ell;
            boundarymat = [boundarymat; boundarydata];
            regprime(:,degcoeff + jj2) = curvepenaltyder(psell,pselljj,tline);
            if jj2 == 1
                centerprime(:,degcoeff + jj2) = center_regder(psell,pselljj,z_0,[1;0],tline);
            else
                centerprime(:,degcoeff + jj2) = center_regder(psell,pselljj,z_0,[0;1],tline);
            end
        end
        Geometry.beta0 = boundarymat;
        fprimenew = direct_sound_soft_splineRK(Var,Geometry); 
        if length(size(fprimenew)) == 2
            fprime = zeros(1,size(fprimenew,1),size(fprimenew,2));
            fprime(1,:,:) = fprimenew;
        else
            fprime = fprimenew;
        end
        gw = 1/errorongivenpoints(uscat,Var.tau) *  sqrt(Var.tau) .* ( uscat_ell - uscat);
        gwshape = reshape(gw.', [],1);

        Fshape = 1/errorongivenpoints(uscat,Var.tau) *  sqrt(Var.tau)*...
            reshape(permute(fprime,[2 1 3]), [nobs * (Var.M+1), degcoeff+2]);

        reg = curvepenalty(psell, tline);
        cen = center_reg(psell,z_0,tline);
        rhs = [-gwshape; -regpar*reg.'; -cenpar*cen];

        Amat = [Fshape.', regpar*regprime.', cenpar*centerprime.'] ...
            * [Fshape; regpar*regprime; cenpar*centerprime];
        bmat = [Fshape.', regpar*regprime.', cenpar*centerprime.'] * rhs;
        update = Amat\bmat;

        updaterad = update(1:end-2);
        updatecen = update(end-1:end);
        errN(j) = norm(gwshape)^2;
        errD(j) = regpar^2 * norm(reg)^2;
        errU(j) = cenpar^2 * norm(cen)^2;
        errG(j) = norm(gwshape)^2 + regpar^2 * norm(reg)^2 + cenpar^2 * norm(cen)^2;


        %% finding the minimum in the direction update using golden ratio search
        smax = 1;
        if golden_section == 1
            golden = 2/(1+sqrt(5));
            X1 = rvecell;
            X1cen = z_0;
            X4 = X1 + smax*updaterad.';
            X4cen = z_0 + smax*updatecen;
            while min(X4)<0
                smax = .9*smax;
                X4 = X1 + smax*updaterad.';
                X4cen = z_0 + smax*updatecen;
            end
            assert(min(X4)>=0)
            X2 = X4 - golden*(X4-X1);
            X2cen = X4cen - golden*(X4cen-X1cen);
            X3 = X1 + golden*(X4-X1);
            X3cen = X1cen + golden*(X4cen-X1cen);
            %prepare func2
            pointsX2 = X2.*[cos(part);sin(part)] + X2cen;
            pointsX2 = [pointsX2, pointsX2(:,1)];
            psX2 = generatesplinehandles(pointsX2);
            Pp = 1/2; Pm = 1/2;
            gX2 = starshapespline(Var.N,0,psX2);
            gpX2 = starshapespline(Var.N,1/6,psX2);
            gmX2 = starshapespline(Var.N,-1/6,psX2);

            uincp = planewaveRK(signal,signalp,dir,Var.Tlag,gpX2.midpt,gpX2.normal,Var.T,Var.M,1,c);
            uincm = planewaveRK(signal,signalp,dir,Var.Tlag,gmX2.midpt,gmX2.normal,Var.T,Var.M,1,c);
            beta0 = (Pp*uincp+Pm*uincm);
            Geometry.ps = psX2;
            Geometry.mindist = minimaldistance(Var,Geometry);
            Geometry.beta0 = beta0;
            uscatX2 = direct_sound_soft_splineRK(Var,Geometry);

            regX2 = curvepenalty(psX2, tline);
            cenX2 = center_reg(psX2,X2cen,tline);
            funct2 = (errorongivenpoints(uscatX2 - uscat,Var.tau) / errorongivenpoints(uscat,Var.tau))^2 + regpar^2 * norm(regX2)^2 + cenpar^2 * norm(cenX2)^2;
            %prepare func3
            pointsX3 = X3.*[cos(part);sin(part)] + X3cen;
            pointsX3 = [pointsX3, pointsX3(:,1)];
            psX3 = generatesplinehandles(pointsX3);
            Pp = 1/2; Pm = 1/2;
            gX3 = starshapespline(Var.N,0,psX3);
            gpX3 = starshapespline(Var.N,1/6,psX3);
            gmX3 = starshapespline(Var.N,-1/6,psX3);

            uincp = planewaveRK(signal,signalp,dir,Var.Tlag,gpX3.midpt,gpX3.normal,Var.T,Var.M,1,c);
            uincm = planewaveRK(signal,signalp,dir,Var.Tlag,gmX3.midpt,gmX3.normal,Var.T,Var.M,1,c);
            beta0 = (Pp*uincp+Pm*uincm);
            Geometry.ps = psX3;
            Geometry.mindist = minimaldistance(Var,Geometry);
            Geometry.beta0 = beta0;
            uscatX3 = direct_sound_soft_splineRK(Var,Geometry);
            regX3 = curvepenalty(psX3, tline);
            cenX3 = center_reg(psX3,X3cen,tline);
            funct3 = (errorongivenpoints(uscatX3 - uscat,Var.tau) / errorongivenpoints(uscat,Var.tau))^2 + regpar^2 * norm(regX3)^2 + cenpar^2 * norm(cenX3)^2;

            %count how many times the farther point is chosen
            counter = 0;
            for kk = 1:numsplits
                if funct2 <= funct3
                    X4 = X3; X4cen = X3cen;
                    X3 = X2; X3cen = X2cen;
                    X2 = X4 - golden * (X4-X1); X2cen = X4cen - golden * (X4cen-X1cen);
                    Xnew = X2; Xnewcen = X2cen;
                    funct3 = funct2;
                    ind = 2;
                else
                    X1 = X2; X1cen = X2cen;
                    X2 = X3; X2cen = X3cen;
                    X3 = X1 + golden*(X4-X1); X3cen = X1cen + golden*(X4cen - X1cen);
                    Xnew = X3; Xnewcen = X3cen;
                    funct2 = funct3;
                    ind = 3;
                    counter = counter +1;
                end
                %prepare funkt
                pointsXnew = Xnew.*[cos(part);sin(part)] + Xnewcen;
                pointsXnew = [pointsXnew, pointsXnew(:,1)];
                psXnew = generatesplinehandles(pointsXnew);
                Pp = 1/2; Pm = 1/2;
                gXnew = starshapespline(Var.N,0,psXnew);
                gpXnew = starshapespline(Var.N,1/6,psXnew);
                gmXnew = starshapespline(Var.N,-1/6,psXnew);

                uincp = planewaveRK(signal,signalp,dir,Var.Tlag,gpXnew.midpt,gpXnew.normal,Var.T,Var.M,1,c);
                uincm = planewaveRK(signal,signalp,dir,Var.Tlag,gmXnew.midpt,gmXnew.normal,Var.T,Var.M,1,c);
                beta0 = (Pp*uincp+Pm*uincm);
                Geometry.ps = psXnew;
                Geometry.mindist = minimaldistance(Var,Geometry);
                Geometry.beta0 = beta0;
                uscatXnew = direct_sound_soft_splineRK(Var,Geometry);
                regXnew = curvepenalty(psXnew, tline);
                cenXnew = center_reg(psXnew,Xnewcen,tline);
                funkt = (errorongivenpoints(uscatXnew - uscat,Var.tau) / errorongivenpoints(uscat,Var.tau))^2 + regpar^2 * norm(regXnew)^2 + cenpar^2 * norm(cenXnew)^2;
                val1 = (errorongivenpoints(uscatXnew - uscat,Var.tau) / errorongivenpoints(uscat,Var.tau))^2;
                val2 = regpar^2 * norm(regXnew)^2;
                val3 = cenpar^2 * norm(cenXnew)^2;
                display([funkt, val1, val2, val3])
                % display([funkt, val1, val2])
                if ind == 2
                    funct2 = funkt;
                else
                    funct3 = funkt;
                end
            end
        else
            X1 = rvecell;
            X1cen = z_0;
            smax = 1;
            X4 = X1 + smax*updaterad.';
            X4cen = z_0 + smax*updatecen;
            while min(X4)<0
                smax = .9*smax;
                X4 = X1 + smax*updaterad.';
                X4cen = z_0 + smax*updatecen;
                fprintf("decrease smax \n")
            end
            assert(min(X4)>=0)
            Xnew = X4;
            Xnewcen = X4cen;
            pointsXnew = Xnew.*[cos(part);sin(part)] + X4cen;
            pointsXnew = [pointsXnew, pointsXnew(:,1)];
            psXnew = generatesplinehandles(pointsXnew);
        end
        total_movement = sum(sqrt(sum((fnval(psell,tline/(2*pi)*psell.breaks(end))-fnval(psXnew,tline/(2*pi)*psell.breaks(end))).^2,1)));
        relative_movement = sum(sqrt(sum((fnval(psell,tline/(2*pi)*psell.breaks(end))-fnval(psXnew,tline/(2*pi)*psXnew.breaks(end))).^2,1)))...
            /sum(sqrt(sum((fnval(psell,tline/(2*pi)*psell.breaks(end))).^2,1)));
        display([relative_movement])

        psell = psXnew;
        rvecell = Xnew;
        z_0 = Xnewcen;
        ss = fnval(psell,tline/(2*pi)*psell.breaks(end));
        ssvec(:,:,j+1) = ss;
        g = starshapespline(Var.N,0,psell);
        gp = starshapespline(Var.N,1/6,psell);
        gm = starshapespline(Var.N,-1/6,psell);

        % 2. step : Solve the forward problem
        % generate boundary data
        Pp = 1/2; Pm = 1/2;
        uincp = planewaveRK(signal,signalp,dir,Var.Tlag,gp.midpt,gp.normal,Var.T,Var.M,1,c);
        uincm = planewaveRK(signal,signalp,dir,Var.Tlag,gm.midpt,gm.normal,Var.T,Var.M,1,c);
        beta0 = (Pp*uincp+Pm*uincm);
        Geometry.ps = psell;
        Geometry.mindist = minimaldistance(Var,Geometry);
        Geometry.beta0 = beta0;
        [uscat_ell,lambdastages] = direct_sound_soft_splineRK(Var,Geometry);
        psder = fnder(psell);
        psderval = fnval(psder, tline/(2*pi)*psder.breaks(end));
        unu_ell =  Var.N/(2*pi)*(1./sqrt(sum(psderval.^2,1))).'.*lambdastages;
        uscatVar{j+1} = uscat_ell;
        j = j+1;
        %     end
    end
    %%
    if save_results == 1
	    clearvars boundarydata boundarymat fprimenew Fshape fprime beta0 lambdastages uincm uincp unu_ell uscat uscat_ell uscatVar
        save(strcat('Paper1/',name,'.mat'));
        clearvars -except rads it
    end
end
%%
plot_pictures = 1;
if plot_pictures
    load('Paper1/Paper11_1.mat')
    PlotIterates(vals,ssvec,Geometry.observationpoints,[1,3,6,9,11,size(ssvec,3)],'Ex1_');
end

