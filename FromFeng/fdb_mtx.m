function bootPxys = fdb_mtx(Pxys_epoch, params)  
% boostrap for the multitaper cross-spectral matrics, resample across tapers
% and time.
%function bootPxys = fdb_mtx(Pxys_epoch, Pxys_t_epoch, params, curf)  

% just compute the bootstrapping for the Pxy. 52-64.
Nfreqs = size(Pxys_epoch,1);
    Ntapers = size(Pxys_epoch,2);
    Nchans = size(Pxys_epoch,3);
    Nt = size(Pxys_epoch,5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Niter =20;

    disp(['*** Using ' num2str(Niter) ' bootstraps'])

    %% bootstrap cross spectra

    disp('- Bootstrapping cross-spectral matrices')

    bootPxys = nan([Nfreqs Nchans Nchans Niter],'like',complex(0));
    for ff = 1:Nfreqs
        if mod(ff,15) == 0
            disp(['    f = ' num2str(params.f(ff))])
        end
        bootPxys(ff,:,:,:) = bootstrap_mtxsptrm(squeeze(Pxys_epoch(ff,:,:,:,:)),Niter);
    end

    clear Pxys_epoch
end

    
    
    
    
    
  