%
function [bootCoh, bootVarexp] = recon_mtx(output, modes)
    %% reconstruct alpha cross-spectrum
if nargin <2
    modes =1:3;
end
Nchans = size(output.eigVals, 1);
Niter = size(output.eigVals, 2);

bootCoh = nan(Nchans,Niter);
bootVarexp = zeros(Niter,1);
for ii =1:Niter
    disp('- Reconstructing single decomposition coherent power')
    disp(['    Using ' num2str(length(modes)) ' modes'])
    U = squeeze(output.eigVects(:,:, ii));
    S = diag(output.eigVals(:, ii));
    [Coh, Varexp] = reconstruct_power(U,S,modes);
    
    bootCoh(:,ii) = Coh;
    bootVarexp(ii) = Varexp;
    
end


end


 