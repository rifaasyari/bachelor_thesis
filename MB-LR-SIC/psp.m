function [P_l] = psp(size, l, L)
% PSP Implements Pre-Stored Patterns (PSP) algorithm proposed in [1]. 
%
%   Input
%       size: Number of columns in the permutation matrix P_l
%       l: Current branch to calculate a permutation matrix for, 2 <= l <=
%           size
%       L: Number of all branches to create patterns for
%
%   Output
%       P_l: Permutation matrix for branch l
%   
%   Examples (here with size = L):
%       P_2 = 0...1     P_3 = 1...0 0   P_L = 1 0 ... 0
%             .  ..           .     1         0 1 ... 0
%             . . .           .   . .         .  .    .
%             .1  .           .  .  .         .    1  .
%             1...0           0 1 . 0         .     0 1
%                                             0 . . 1 0  
%
%   [1]: R. Fa and R.C. de Lamare, "Design of Adaptive Multi-Branch SIC
%       Receivers for MIMO Spatial Multiplexing Systems", in International
%       Symposium on Wireless Communications Systems (ISWCS),
%       Siena-Tuscany, Italy, 2009, pp.575 - 579
%
%   See also CNBO
    
    s = floor((l - 2) * size / L);
    P_l = [eye(s) zeros(s,size - s); zeros(size - s,s) rot90(eye(size - s))];
    
    addflops(2 + flops_div() + 3);
end
