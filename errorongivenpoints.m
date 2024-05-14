function err = errorongivenpoints(data,tau)
%%%%%%%%%%%%%%
%% Input
%  data: the vector whose norm should be returned
%  tau: the scaling parameter T/M!
%% Output
%  err : the norm of data.
%%%%%%%%%%%%%%

err = sqrt(tau*sum(sum(abs(data).^2)));

end