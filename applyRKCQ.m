function [u,eta] = applyRKCQ(F,g,d2,T,N,A,b,c)
%%%%%%%%%%%%%%%%%%%%%
% Input: 
% F: time-harmonic function handle F(s,b): C x C^d1 -> C^d2, d1 number of
%    spatial points
% g: real-valued d1 x (m*N + 1)-matrix, (m amount of stages of RK)
% d2: number of observation points
% T: final time
% N: amount of timesteps
% A,b,c: coefficients from the RKCQ
%
% Output: 
% u: the scattered field at the observation points. Depending on the input,
%    this may also be the time-dependent domain derivative u' 
% eta: this is the density satisfying V*eta = beta0. eta is required for
%      the evaluation of the domain derivative, since it may approximate
%      \partial_\nu (u+u^i) up to some scaling.
%%%%%%%%%%%%%%%%%%%%%%

dt = T/N;
m = length(c);
tol = 10^(-14);
rho = tol^(1/(2*N));
L = N;
d1 = length(g(:,1));
g_fft = zeros(d1,m*L);
phi_inter = zeros(d1,L);
if nargout == 2
    %% Step 1: Preprocessing
    for stageInd=1:m
        phi_inter(:,1:N) = rho.^(linspace(0,N-1,N)).*g(:,stageInd+1:m:m*N+1);
        g_fft(:,stageInd:m:m*L) = fft(phi_inter,L,2);
    end
    %% Step 2: The time-harmonic operators
    s_vect = rho*exp(-1j*2*pi*linspace(0,L-1,L)/L);
    HalfL = ceil(L/2);
    % phi_hat   = zeros(d2,m*L);
    phi_hatVec = zeros(d2,m,HalfL+1);
    eta_hatVec = zeros(size(g,1),m,HalfL+1);
    % avoid unnecessary communication overhead! split the broadcast
    % variable g_fft before!
    if m ==2
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end);
            G = [g_fft1;g_fft2];
        elseif m == 3
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end); g_fft3 = g_fft(:, 3:m:end);
            G = [g_fft1;g_fft2;g_fft3];
        elseif m ==4
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end); g_fft3 = g_fft(:, 3:m:end); g_fft4 = g_fft(:, 4:m:end);
            G = [g_fft1;g_fft2;g_fft3;g_fft4];
    end
    % vectorized
    parfor j = 1:HalfL+1
        deltaMatrix = inv(A + s_vect(j) * 1.0 / (1 - s_vect(j)) * ones(m, 1) * b');
        [T, deltaEigsM] = eig(deltaMatrix / dt);
        deltaEigs = diag(deltaEigsM);
        Tinv = inv(T);
        rhsStages = reshape(G(:,j),d1,m)*(Tinv.');
        lhsStagesVec = [];
        etaVec = [];
        for stageInd =1:m
            [SVinvrhs,Vinvrhs] = F(deltaEigs(stageInd),rhsStages(:,stageInd));
            lhsStagesVec = [lhsStagesVec,SVinvrhs];
            etaVec = [etaVec, Vinvrhs];
        end
        phi_hatVec(:,:,j) = lhsStagesVec * T.';
        eta_hatVec(:,:,j) = etaVec * T.';
    end
    phi_hat = reshape(phi_hatVec,d2,(HalfL+1)*m);
    eta_hat = reshape(eta_hatVec,size(g,1),(HalfL+1)*m);
    %%
    % Mirroring the second part of the frequencies by complex conjugation (vectorized)
    freqInd = (HalfL+1):(L-1);
    sourceIndices = bsxfun(@plus,m * (L - freqInd) , (1:m).');
    targetIndices = bsxfun(@plus,freqInd  * m , (1:m)');
    % Perform vectorized complex conjugation and assignment
    phi_hat(:, targetIndices) = conj(phi_hat(:, sourceIndices));
    eta_hat(:, targetIndices) = conj(eta_hat(:, sourceIndices));
    %% Step 3: Postprocessing
    %reminder:
    %d1 = length(g(:,1)); phi_inter = zeros(d2,L);
    phi_solL = zeros(d2,m*L);

    for stageInd=1:m
        phi_inter = real(ifft(phi_hat(:,stageInd:m:m*L),L,2));
        phi_solL(:,stageInd:m:m*L) = rho.^(-linspace(0,L-1,L)).*phi_inter;

        eta_inter = real(ifft(eta_hat(:,stageInd:m:m*L),L,2));
        eta_solL(:,stageInd:m:m*L) = rho.^(-linspace(0,L-1,L)).*eta_inter;
    end
    phi_sol = zeros(d2,m*N+1);
    phi_sol(:,2:m*N+1) = phi_solL(:,1:m*N);

    u=phi_sol;

    eta_sol = zeros(size(g,1),m*N+1);
    eta_sol(:,2:m*N+1) = eta_solL(:,1:m*N);

    eta=eta_sol;
else
    spacepoints = F(1,1);
    if spacepoints == d1 % this is the usual situation.
        %% Step 1: Preprocessing
        for stageInd=1:m
            phi_inter(:,1:N) = rho.^(linspace(0,N-1,N)).*g(:,stageInd+1:m:m*N+1);
            g_fft(:,stageInd:m:m*L) = fft(phi_inter,L,2); % it's the same as fft(phi_inter,[],2) 200x400
        end
        %% Step 2: The time-harmonic operators
        s_vect = rho*exp(-1j*2*pi*linspace(0,L-1,L)/L);
        HalfL = ceil(L/2);
        phi_hatVec = zeros(d2,m,HalfL+1);
        % avoid unnecessary communication overhead! split the broadcast
        % variable g_fft before!
        if m ==2
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end);
            G = [g_fft1;g_fft2];
        elseif m == 3
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end); g_fft3 = g_fft(:, 3:m:end);
            G = [g_fft1;g_fft2;g_fft3];
        elseif m ==4
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end); g_fft3 = g_fft(:, 3:m:end); g_fft4 = g_fft(:, 4:m:end);
            G = [g_fft1;g_fft2;g_fft3;g_fft4];
        end
        % vectorized
        parfor j = 1:HalfL+1
            deltaMatrix = inv(A + s_vect(j) * 1.0 / (1 - s_vect(j)) * ones(m, 1) * b');
            [T, deltaEigsM] = eig(deltaMatrix / dt);
            deltaEigs = diag(deltaEigsM);
            Tinv = inv(T);
            rhsStages = reshape(G(:,j),d1,m)*(Tinv.');
            lhsStagesVec = [];
            for stageInd =1:m
                SVinvrhs = F(deltaEigs(stageInd),rhsStages(:,stageInd));
                lhsStagesVec = [lhsStagesVec,SVinvrhs];
            end
            phi_hatVec(:,:,j) = lhsStagesVec * T.';
        end
        phi_hat = reshape(phi_hatVec,d2,(HalfL+1)*m);
        %%
        % Mirroring the second part of the frequencies by complex conjugation (vectorized)
        freqInd = (HalfL+1):(L-1);
        sourceIndices = bsxfun(@plus,m * (L - freqInd) , (1:m).');
        targetIndices = bsxfun(@plus,freqInd  * m , (1:m)');
        % Perform vectorized complex conjugation and assignment
        phi_hat(:, targetIndices) = conj(phi_hat(:, sourceIndices));
        %% Step 3: Postprocessing
        %reminder: 
        %d1 = length(g(:,1)); phi_inter = zeros(d2,L);
        phi_solL = zeros(d2,m*L);

        for stageInd=1:m
            phi_inter = real(ifft(phi_hat(:,stageInd:m:m*L),L,2));
            phi_solL(:,stageInd:m:m*L) = rho.^(-linspace(0,L-1,L)).*phi_inter;
        end
        phi_sol = zeros(d2,m*N+1);
        phi_sol(:,2:m*N+1) = phi_solL(:,1:m*N);

        u=phi_sol;

        eta=[];
    else % this is only for the efficient evaluation of the Frechet derivative.
        % note that we perform the very same steps as before but with 
        % Q*d2 instead of d2. In the end we reshape in order to get a 
        % d2,m*N+1,Q matrix. So, all the corr. directions are found in the
        % third dimension.
        Q = d1/spacepoints;
        %% Step 1: Preprocessing
        for stageInd=1:m
            phi_inter(:,1:N) = rho.^(linspace(0,N-1,N)).*g(:,stageInd+1:m:m*N+1);
            g_fft(:,stageInd:m:m*L) = fft(phi_inter,L,2); %% Hier anfangen!!!! ok.
        end
        %% Step 2: The time-harmonic operators
        s_vect = rho*exp(-1j*2*pi*linspace(0,L-1,L)/L);
        HalfL = ceil(L/2);
        u = zeros(Q*d2,m,HalfL+1);
        % avoid unnecessary communication overhead! split the broadcast
        % variable g_fft before!
        if m ==2
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end);
            G = [g_fft1;g_fft2];
        elseif m == 3
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end); g_fft3 = g_fft(:, 3:m:end);
            G = [g_fft1;g_fft2;g_fft3];
        elseif m ==4
            g_fft1 = g_fft(:, 1:m:end); g_fft2 = g_fft(:, 2:m:end); g_fft3 = g_fft(:, 3:m:end); g_fft4 = g_fft(:, 4:m:end);
            G = [g_fft1;g_fft2;g_fft3;g_fft4];
        end
        % vectorized
        parfor j = 1:HalfL+1
            deltaMatrix = inv(A + s_vect(j) * 1.0 / (1 - s_vect(j)) * ones(m, 1) * b');
            [T, deltaEigsM] = eig(deltaMatrix / dt);
            deltaEigs = diag(deltaEigsM);
            Tinv = inv(T);
            rhsStages = reshape(G(:,j),d1,m)*(Tinv.');
            lhsStagesVec = [];
            for stageInd =1:m
                rhsStagesShape = reshape(rhsStages(:,stageInd),spacepoints,Q); 
                SVinvrhsShape = F(deltaEigs(stageInd),rhsStagesShape);
                SVinvrhs = reshape(SVinvrhsShape,Q*d2,1);
                lhsStagesVec = [lhsStagesVec,SVinvrhs];
            end
            u(:,:,j) = lhsStagesVec * T.';
        end
        %% this is replaced by the vectorization below.
        % phi_hat = [];
        % for kk = 1 : Q
        %     phi_hat = [phi_hat; reshape(u((kk-1)*d2+1:kk*d2,:,:),d2,(HalfL+1)*m)];
        % end
        %%
        % Calculate the linear indices for extraction
        indices = reshape(((1:Q)' - 1) * d2 + 1, 1, 1, Q) + (0:d2-1);
        % Extract the blocks and reshape 
        u = reshape(u(indices, :, :), Q * d2, (HalfL + 1) * m);
        % to save space: u = phi_hat at this point.
        %%
        % Mirroring the second part of the frequencies by complex conjugation (vectorized)
        freqInd = (HalfL+1):(L-1);
        sourceIndices = bsxfun(@plus,m * (L - freqInd) , (1:m).');
        targetIndices = bsxfun(@plus,freqInd  * m , (1:m)');
        % Perform vectorized complex conjugation and assignment
        u(:, targetIndices) = conj(u(:, sourceIndices));
        %% Step 3: Postprocessing
        % Vectorization. We overwrite phi_hat_Vec all the time to
        % save memory
        stageIndices = bsxfun(@plus, 1:m:m*L , (0:m-1).');
        % I want to do a block wise reshapement here. This does the trick.
        u = reshape(reshape(u(:,stageIndices),Q*d2,L,m),[],L); % this works
        u = real(ifft(u, [], 2));
        u = rho.^(-linspace(0, L-1, L)) .* u;%bsxfun(@times, rho.^(-linspace(0, L-1, L)), phi_inter2);
        % reshape the block again
        u = reshape(reshape(u,Q*d2,m,L),[],L*m);
        %%        
        % this here puts zeros to the first column of u. This corresponds
        % to the value at time 0, that is not treated throughout the RKCQ!
        u(:,end+1) = 0;
        u = circshift(u,1,2);
        % Reshape phi_sol into a 3D matrix
        u = reshape(u, d2, Q, m*N+1);

        % Use indexing to assign values directly
        % u(:, :, 1:Q) = permute(u, [1, 3, 2]);
        u = permute(u, [1, 3, 2]);
        eta=[];
    end
end
end