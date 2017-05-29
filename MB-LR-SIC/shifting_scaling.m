function [estimate] = shifting_scaling(z, n, kappa, T_inv)
% SHIFTING_SCALING Shifting and scaling operations to solve quantisation
% problem of lattice-reduction aided MIMO detection [1].
%
%   Input
%       z: Detected signal in the reduced lattice space
%       n: Index of the detected component
%       kappa: Complex offset in the normalised constellation space. 
%           - kappa = 1+j for all M-QAM constellations
%           - kappa = 1 for all PAM constellations 
%       T_inv: Inverse of the uni-modular matrix used for lattice reduction
%
%   Output
%       estimate: Estimation of transmitted signal ẑ in the reduced
%           lattice. Ensures that ẑ stems from lattice reduced
%           constellation with consecutive integer symbols. 
%
%   [1]: D. Milford and M. Sandell, "Simplified Quantisation in a
%       Reduced-Lattice MIMO Decoder", IEEE Communication Letters, vol. 15,
%       no. 7, pp. 725 - 727, July 2011
%
%   See also MR_LR_SIC

    cmpl_lattice_offset = kappa * T_inv(n,:) * ones(length(T_inv(n,:)), 1);
    estimate = 2 * round((z(1) - cmpl_lattice_offset)/2) + ... 
        cmpl_lattice_offset;
    
    addflops(2 * flops_mul(1, length(T_inv(n,:)), 1) + 5 + ... 
        flops_div() + 2);

end
