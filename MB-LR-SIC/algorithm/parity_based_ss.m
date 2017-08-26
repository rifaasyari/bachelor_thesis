function [estimate] = parity_based_ss(z, n, kappa, T_inv, var_s, modulation)
% PARITY_BASED_SS Alternative procedure to shifting and scaling operations.
% The algorithm is taken from [1].
%
%   Input
%       z (scalar): Detected signal component in the reduced lattice space
%       n (scalar): Index of the detected component
%       kappa (scalar): Complex offset in the normalised constellation space. 
%           - kappa = 1+j for all M-QAM constellations
%           - kappa = 1 for all PAM constellations 
%       T_inv (matrix): Inverse of the uni-modular matrix used for lattice 
%           reduction
%       var_s (scalar): Symbol variance
%       modulation (string): Used modulation: "QPSK" or "16-QAM"
%
%   Output
%       estimate (scalar): Estimation of transmitted signal s in the reduced
%           lattice. Ensures that s stems from lattice reduced
%           constellation with consecutive integer symbols. 
%
%   [1]: D. Milford and M. Sandell, "Simplified Quantisation in a
%       Reduced-Lattice MIMO Decoder", IEEE Communication Letters, vol. 15,
%       no. 7, pp. 725 - 727, July 2011
%
%   See also MR_LR_SIC

    rho = sum(T_inv(n,:));

    if strcmp(modulation, 'QPSK')
        norm = 1/sqrt(var_s);
    elseif strcmp(modulation, '16-QAM')
        d = sqrt(2/5 * var_s);
        norm = 2/d;
    else
        error('Invalid modulation: %s', modulation);
    end
    
    z = z * norm;  % Normalize lattice for Shifting and Scaling
    parity = @(x) (mod(real(x), 2) + 1j*mod(imag(x),2));
    cxor = @(x,y) (xor(real(x),real(y)) + 1j*xor(imag(x),imag(y)));  
    % complex XOR
    
    mulmod2 = @(x,y) (xor(and(real(x),real(y)), and(imag(x), imag(y))) + ...
        1j*xor(and(real(x),imag(y)), and(imag(x),real(y))));
        
    estimate = floor(z) + cxor(parity(floor(z)), mulmod2(kappa, parity(rho)));
    
    addflops(flops_mul(1, length(T_inv(n,:)), 1) + 5 + ... 
        flops_div() + 2);
    
    estimate = estimate / norm;  % Scale back 
    
end
