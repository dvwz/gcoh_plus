function output = eig_gcoh(Pxys_epoch,params, curf) 
  
% Pxys_epoch Nfreq x Nchans x Nachans x Niter

Niter = size(Pxys_epoch, 4);
Nchans = size(Pxys_epoch, 2);
    [~,fpos] = min(abs(params.f - curf));
    findx = params.findx(fpos);
    assert(length(findx) == 1)

    Pxys_alpha = squeeze(Pxys_epoch(findx,:,:,:,:)); % now Nchans x Nchans x Niter

    % mean over tapers
    %medPxys = Pxys_alpha; %squeeze(mean(Pxys_alpha,1)); % now Nchans x Nchans x Nchunks
    % nanmedian over time, which means nans are treated as missing
    %medPxys = squeeze(nanmedian(real(medPxys),3))+...
    %    1i*squeeze(nanmedian(imag(medPxys),3));
    eigVects =  zeros(Nchans, Nchans, Niter);
    eigVals  =  zeros(Nchans, Niter);
    for ii = 1:Niter
    %compute eig decomposition
        [U,S]=eignan_gcoh(squeeze(Pxys_alpha(:,:,ii)));
        [~, sv_ind]=sort(diag(S),'descend'); %we will plot the dominant singval
        S = S(sv_ind,sv_ind);
        U = U(:,sv_ind); % put eigenvectors in appropriate order
        
        eigVects(:,:, ii) = U;
        
        eigVals(:, ii)  = diag(S);
        
    end
    
    
    output.(['eigVects']) = eigVects;
    output.(['eigVals']) = eigVals;

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% BEGIN BOOTSTRAP RUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     
% %% bootstrap eignan
% 
% disp('- Bootstrapping eigendecompositions')
% 
% bootU = nan([Nfreqs Nchans Nchans Niter]);
% bootS = nan([Nfreqs Nchans Niter]);
% for ff = 1:Nfreqs
%     %disp(['    f = ' num2str(params.f(ff))])
%     [bU, bS] = bootstrap_eignan(squeeze(bootPxys(ff,:,:,:)));
%     bootU(ff,:,:,:) = bU;
%     bootS(ff,:,:) = bS;
% end
% 
% %% bootstrap cohpower reconstructions 
% % throw away nondominant eigvect/values
% 
% disp('- Bootstrapping coherent power')
% 
% bootCoh = nan([Nfreqs Nchans Niter]);
% %bootVarexp = nan([Nfreqs Niter]);
% %bootProj = nan([Nfreqs Nchans Nchans Niter]);
% for ff = 1:Nfreqs
%     %disp(['    f = ' num2str(params.f(ff))])
%     [bCoh, ~, ~] = bootstrap_cohpower(squeeze(bootU(ff,:,:,:)), squeeze(bootS(ff,:,:)), modes);
%     bootCoh(ff,:,:) = bCoh;
%     %bootVarexp(ff,:) = bVarexp;
%     %bootProj(ff,:,:,:) = bProj;
% end
% 
% %% assign outputs
% %output.(['curepoch']) = curepoch;
% 
% output.(['params']) = params;
% 
% output.(['t_for_Pxys']) = Pxys_t_epoch;
% output.(['findx']) = findx;
% 
% output.(['medPxys'])=medPxys;
% output.(['singvals_sorted'])=sv_srtd;
% %output.(['Pxys']) = allPxys;
% output.(['U'])=U;
% output.(['S'])=S;
% output.(['Coh']) = Coh;
% %output.(['Varexp']) = Varexp;
% %output.(['Proj']) = Proj;
% 
% %output.(['bootPxys']) = bootPxys;
% output.(['bootU']) = bootU;
% output.(['bootS']) = bootS;
% output.(['bootCoh']) = bootCoh;
% %output.(['bootVarexp']) = bootVarexp;
% %output.(['bootProj']) = bootProj;
% output.(['niter']) = Niter;
% 
% %output.(['gcoh_labels'])=gcoh_labels;
% 
% output.(['curf'])=curf;
% %output.(['loctime']) = loctime;
% 
% %output.(['artdir']) = artdir;