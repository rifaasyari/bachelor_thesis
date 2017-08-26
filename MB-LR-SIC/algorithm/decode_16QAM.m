function [b] = decode_16QAM(s, mod)
% DECODE_16QAM Decode 16-QAM symbols into four bits for each symbol. Symbols
% are supposed to be Gray coded in the following manner
%        1       2        3       4
%
%    1   *       *    ^   *       *
%       1011    1010  |  0010    0011
%                     |
%    2   *       *    |   *       *
%       1001    1000  |  0000    0001
%      -------------(0,0)------------->
%    3   *       *    |   *       *
%       1101    1100  |  0100    0101
%                     |
%    4   *       *    |   *       *
%       1111    1110  |  0110    0111
% 
% INPUT s (scalar or vector): Sequence of 16-QAM symbols
%       mod (matrix): 16-QAM symbol constellation
% 
% OUTPUT b (scalar or vector): Bit sequence of symbols s of four times length
%
% See also DECODE_QPSK

    b = zeros(1, 4*length(s));
    for k = 1:length(s)
        pos = (k-1) * 4 + 1;
        S = repmat(s(k), size(mod));
        ml = abs(mod-S).^2;
        [~, min_rows] = min(ml);
        [~, min_col] = min(min(ml));
        min_row = min_rows(min_col);
        switch min_row
            case 1
                switch min_col
                    case 1
                        b(pos:pos+3) = [1 0 1 1];
                    case 2
                        b(pos:pos+3) = [1 0 1 0];
                    case 3
                        b(pos:pos+3) = [0 0 1 0];
                    case 4
                        b(pos:pos+3) = [0 0 1 1];
                end
            case 2
                switch min_col
                    case 1
                        b(pos:pos+3) = [1 0 0 1];
                    case 2
                        b(pos:pos+3) = [1 0 0 0];
                    case 3
                        b(pos:pos+3) = [0 0 0 0];
                    case 4
                        b(pos:pos+3) = [0 0 0 1];
                end
            case 3
                switch min_col
                    case 1
                        b(pos:pos+3) = [1 1 0 1];
                    case 2
                        b(pos:pos+3) = [1 1 0 0];
                    case 3
                        b(pos:pos+3) = [0 1 0 0];
                    case 4
                        b(pos:pos+3) = [0 1 0 1];
                end
            case 4
                switch min_col
                    case 1
                        b(pos:pos+3) = [1 1 1 1];
                    case 2
                        b(pos:pos+3) = [1 1 1 0];
                    case 3
                        b(pos:pos+3) = [0 1 1 0];
                    case 4
                        b(pos:pos+3) = [0 1 1 1];
                end
        end
    end

end
