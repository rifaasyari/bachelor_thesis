function [b] = decode_qpsk(s, ~)
% DECODE_QPSK Decode QPSK symbols into two bits for each symbol. Symbols
% are supposed to be Gray coded in the following manner
%   1+1j -> '00', -1+1j -> '01', -1-1j -> '11', 1-1j -> '10'
% 
% INPUT s (scalar or vector): Sequence of QPSK symbols
% 
% OUTPUT b (scalar or vector): Bit sequence of symbols s of double length

    b = zeros(1, 2*length(s));
    for k = 1:length(s)
        if abs(s(k)) == 0  % s(k) = 0 + 0j
            b(k:k+1) = randi([0,1], 1, 2);  % Random bits
            % fprintf('Guessing random bits...\n');
        else
            if real(s(k)) >= 0 && imag(s(k)) >= 0
                b(k:k+1) = [0 0];
            elseif real(s(k)) <= 0 && imag(s(k)) >= 0
                b(k:k+1) = [0 1];
            elseif real(s(k)) <= 0 && imag(s(k)) <= 0
                b(k:k+1) = [1 1];
            elseif real(s(k)) >= 0 && imag(s(k)) <= 0
                b(k:k+1) = [1 0];
            else
                disp(s(k));
                error('Invalid symbol');
            end
        end
    end

end
