function [z_hat] = shifting_scaling(z, kappa, T_inv, n, var_s, ... 
    modulation)
% SHIFTING_SCALING Shifting and scaling operations to solve quantisation
% problem of lattice-reduction aided MIMO detection [1].
%
%   Input
%       z (scalar): Detected signal component in the reduced lattice space
%       kappa (scalar): Complex offset in the normalised constellation space. 
%           - kappa = 1+j for all M-QAM constellations
%           - kappa = 1 for all PAM constellations 
%       T_inv (matrix): Inverse of uni-modular transformation matrix used for 
%           lattice reduction
%       n (scalar): Index of detected signal component
%       var_s (scalar): Symbol variance
%       modulation (string): Used modulation: "QPSK" or "16-QAM"
%       normalize(bool): If true, normalize symbols. Otherwise not. 
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
    
    cmpl_lattice_offset = kappa * sum(T_inv(n,:));

    if strcmp(modulation, 'QPSK')
        a = sqrt(2 * var_s);
        norm_factor = 2/a;
    elseif strcmp(modulation, '16-QAM')
        d = sqrt(2/5 * var_s);
        norm_factor = 2/d;
    else
        error('Invalid modulation: %s', modulation);
    end
    
    z = z * norm_factor;  % Normalize lattice for Shifting and Scaling

    z_hat = 2 * round((z - cmpl_lattice_offset)/2) + cmpl_lattice_offset;

        
%     addflops(flops_mul(1, length(T_inv(n,:)), 1) + 5 + ... 
%         flops_div() + 2);
    
    z_hat = z_hat / norm_factor;  % Scale back
    
end
