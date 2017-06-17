function [mod] = create_16_QAM(var_s)
% CREATE_16_QAM Creates a 16-QAM constellation with given separation
%   between the symbols
%
% INPUT var_s: Symbol variance. For 16-QAM it follows that the symbol
%   separation d is
%   var_s = 5/2 * d^2 <=> d = sqrt(2/5 * var_s);
%
% OUTPUT mod: 16 symbols with separation d

    d = sqrt(2/5 * var_s);
    re = -3/2*d:d:3/2*d;
    im = 1j * (3/2*d:-d:-3/2*d);
    [Re, Im] = meshgrid(re,im);
    mod = Re+Im;
    
    % plot mod
    % re = reshape(Re, [1,16]);
    % im = reshape(-Im*1j, [1,16]);
    % scatter(re, im, 100, 'filled');
end
