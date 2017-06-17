function [mod] = create_MQAM(M, var_s)
% CREATE_MQAM Creates a M-QAM constellation with given separation
%   between the symbols
%
% INPUT M (scalar): Modulation level (e.g. 16, 64, 256). M must be a power 
%           of 2
%       var_s (scalar): Symbol variance. For M-QAM it follows that the 
%           symbol separation d is
%               var_s = 1/6 * (M-1) * d^2 <=> d = sqrt(6 * var_s / (M-1));
%
% OUTPUT mod (matrix): M symbols with separation d

    d = sqrt(6 * var_s / (M-1));
    re = -(sqrt(M)-1)/2*d:d:(sqrt(M)-1)/2*d;
    im = 1j * ((sqrt(M)-1)/2*d:-d:-(sqrt(M)-1)/2*d);
    [Re, Im] = meshgrid(re,im);
    mod = Re+Im;
    
    %% plot mod
    % re = reshape(Re, [1,M]);
    % im = reshape(-Im*1j, [1,M]);
    % scatter(re, im, 50, 'filled');
end
