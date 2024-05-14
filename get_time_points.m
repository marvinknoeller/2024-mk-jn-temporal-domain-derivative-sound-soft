function ts = get_time_points(T,tau, m, c)
% get the time points
%% Input
%  T: final time
%  tau: time step size tau = T/M
%  m: number of RK stages
%  c: coefficient c from RK
%% Output
%  ts: the time points needed for the RKCQ
%%%%%%%%%%%%%%
    N = round(T / tau);
    ts = zeros(1, m * N + 1);
    for j = 0:N-1
        for k = 1:m
            ts(j * m + k + 1) = (j + c(k));
        end
    end
    ts = tau * ts;
end
