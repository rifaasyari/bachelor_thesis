function [H, D, P] = ls_channel_estimation_mimo(y_i, s_i, lambda, delta, D, P)
% LS_CHANNEL_ESTIMATION_MIMO Implements the least square (LS) channel
% estimation algorithm for MIMO systems proposed in [1]. 
% This method assumes a block-flat fading channel with T_B >> T_S where T_B
% and T_S denotes, respectively, the time where the channel is constant and 
% the symbol time. Channel estimation is done at the start of each block by
% sending a sequence of training data to the receiver. The estimation of 
% the channel gain matrix will be used for the remaining data that follows. 
% After T_B time has passed the channel gains will change to uncorrelated
% new values so that the new channel matrix has to be estimated again.
%
%   Input
%       y_i (vector): Received signal for the channel estimation at time
%           instant i of training
%       s_i (vector): Training signal at time instant i of training sent by
%           the transmitter. Training signals must be known to both
%           transmitter and receiver
%       lambda (scalar): Forgetting factor in channel estimation. Typically, 
%           0.5 < lambda < 1. Refer to [1] for a optimum value estimation
%           of lambda.
%       delta (scalar): Small constant used for initialization
%       D (matrix): Matrix to calculate channel estimation H. Must be given
%           from last iteration. If this is the first iteration then it is 
%           initialized by delta
%       P (matrix): Matrix to calculate channel estimation H. Must be given
%           from last iteration. If this is the first iteration then it is
%           initialized as zero matrix.
%
%   Output:
%       H: Channel estimation of i-th iteration
%       D: Matrix D from i-th iteration
%       P: Matrix P from i-th iteration
%
%   [1]: E. Karami, "Tracking Performance of Least Squares MIMO Channel
%       Estimation Algorithm", IEEE Transactions on Communications, vol.
%       55, no. 11, pp. 2201 - 2209, Nov. 2007

    if ~exists('D', 'var') && ~exists('P', 'var')
        N_r = length(y_i);  % Number of total receive antennas
        N_t = length(s_i);  % Number of total transmit antennas 
        D = zeros(N_r, N_t);
        P = 1/delta * eye(N_t);
    end
    
    D = lambda * D + y_i * ctranspose(s_i);
    P = lambda^-1 * P - (lambda^-2 * P * s_i * ctranspose(s_i) * P / ... 
        1 + lambda^-1 * ctranspose(s_i) * P * s_i);
    H = D * P;
end
