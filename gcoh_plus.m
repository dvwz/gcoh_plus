%GCOH_PLUS Top-level function of the gcoh+ toolbox
% Description:
%    Could be a wrapper function, or object-oriented class method
%    Desired functionality: creates output data path, calls all versions of
%    the toolbox (per-window/per-epoch, single/bootstrap/null mtxspctrms,
%    parallel/nonparallel,...)
%    Settings: # of components to compute, multitaper settings
%    
% Usage:
%    [singvals, singvecs_l, singvecs_r, gcohmat, times, freqs, gcohsig, mc_samples] = gcoh_plus(data, Fs, T, mw, mtparams, Niter)
%    [output1,output2] = gcoh_plus(input1,input2)
%
% Inputs:
%    data - time series (in form nsamples x nchannels)
%    epoch - times of epoch segments, in seconds (in the form 
%        [tstart tend], or a nsegs x 2 array, or [])
%    mtparams - structure with fields of multitaper params such as
%        tapers, Fs, fpass, etc
%    gcohparams (optional) - structure with fields savedir, outputdir,
%        version ('epoch', 'window'), inference ('none', 'fdb', 'null'), 
%        parallel ('true' | 'false'), Niter
%
% Outputs:
%    eigenvectors_l - left eigenvectors (in form nchannels x nchannels)
%    eigenvectors_r - right eigenvectors (in form nchannels x nchannels)
%    eigenvalues - eigenvalues (in form nchannels x 1)
%    times - time axis vector for plotting
%    freqs - frequency axis vector for plotting
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Other requirements: none
%
% See also: gcoh.m
%
% Copyright Apr-2020, David Zhou, dwzhou@mit.edu
% Last revision 06-Apr-2020
%------------------------------------------------

function gcohout = gcoh_plus(data, ...
                             epochs, ...
                             mtparams, ...
                             gcohparams)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROCESS INPUTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% MTPARAMS
if ~isstruct(mtparams)
    error('gcoh2:paramsstruct','Error:  params must be a structure.')
elseif ~isfield(mtparams, 'Fs')
    error('gcoh2:missingFs','Error:  Please enter a value for params.Fs')
end

% default multitaper parameters for Chronux
if ~isfield(mtparams, 'pad'); mtparams.pad = 0; end
if ~isfield(mtparams, 'tapers'); mtparams.tapers = [3 2]; end
if ~isfield(mtparams, 'fpass'); mtparams.fpass = [0,(Fs/2)]; end
if ~isfield(mtparams, 'err'); mtparams.err = 0; end
if ~isfield(mtparams, 'trialave'); mtparams.trialave = 0; end
if ~isfield(mtparams, 'Fs');
    mtparams.Fs = input('Fs is not present in mtparams. Please enter sampling rate: ');
end
if ~isfield(mtparams, 'movingwin')
    if Fs <= 25
        mtparams.movingwin = [20 20];
    else
        mtparams.movingwin = [3 3];
    end
end

% GCOHPARAMS
if ~isstruct(gcohparams)
    error('gcoh2:gcohparamsstruct','Error:  gcohparams must be a structure.')
end

if ~isfield(gcohparams, 'savedir'); gcohparams.savedir = pwd; end
if ~isfield(gcohparams, 'outputdir'); gcohparams.outputdir = pwd; end
if ~isfield(gcohparams, 'inference'); gcohparams.inference = 'none'; end
if ~isfield(gcohparams, 'version'); gcohparams.version = 'window'; end
if ~isfield(gcohparams, 'Niter')
    switch gcohparams.inference
        case 'none' % just one svd per mtx
            gcohparams.Niter = 1;
        case 'fdb' % frequency domain bootstrap
            gcohparams.Niter = 200;
        case 'null'
            gcohparams.Niter = 200;
    end
end
if ~isfield(gcohparams, 'parallel'); gcohparams.parallel = 0; end

% check data
mwFs = mtparams.movingwin*mtparams.Fs; % [n samples per window, n samples step size]
assert(size(data,1) > mwFs(1)) % the data almost always has to be larger in samples than n_channels
Nchans = size(data,2); % number of channels

% check T
T_samples = 1:size(data,1); % vector of sample indices
T_sec = T_samples ./ mtparams.Fs; % time vector of samples in seconds
if isempty(epochs)
    T_winstart = T_samples(1:mwFs(2):end-mwFs(1)); % vector of window start indices
    T_secout = T_sec(1:mwFs(2):end-mwFs(1)); % time vector of gcoh outputs in seconds
elseif size(epochs,2) == 2
    T_winstart = [];
    T_secout = [];
    for e = 1:size(epochs,1)
        s_start = ceil(epochs(e,1)*mtparams.Fs); % the sample that the epoch should be started
        s_end = floor(epochs(e,2)*mtparams.Fs); % the sample the the epoch should be ended
        T_winstart = [T_winstart (s_start:mwFs(2):s_end-mwFs(1))];
        T_secout = [T_secout ((s_start:mwFs(2):s_end-mwFs(1)) ./ mtparams.Fs)];
    end
end

Nmw = length(T_winstart); % number of moving windows
times = T_secout;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INIT OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gcohout = struct;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------
% CUT OUT DATA WINDOWS FOR EPOCHS
%---------------------------------------

data_windowed = nan(Nmw*mwFs(1),Nchans);
di = 1;
for e = 1:length(T_winstart)
    data_windowed(di:di+mwFs(1)-1,:) = data(T_winstart(e):T_winstart(e)+mwFs(1)-1,:);
    di = di + mwFs(1);
end
data = data_windowed; clear data_windowed;

%---------------------------------------
% COMPUTE CROSS-SPECTRUM WITH ALL TAPERS
%---------------------------------------

% first chunk of data, linear detrend
data_chunk = detrend(data(1:mwFs(1),:));

fprintf('Computing multitaper cross-spectrum... ');
[Pxy,freqs] = pkmtxsptrm(data_chunk,mtparams);
Nfreqs = numel(freqs); frequencies = 1:Nfreqs;
% preallocate Pxy with Nfreqs x Ntapers x Nchans x Nchans x Nmw
Pxy(end,end,end,end,Nmw) = zeros(1,1,'like',complex(0)); 

% multi-taper spectrum, to infinity and beyond!
for qq = 1:Nmw
    
    qp = qq-1;
    data_chunk = detrend(data(mwFs(2)*qp+(1:mwFs(1)),:));
    Pxy(:,:,:,:,qq) = pkmtxsptrm(data_chunk,mtparams);
    
end
clear data_chunk data;
disp('Done.')

%---------------------------------------
% PRINT WHAT VERSION IS BEING RUN
%---------------------------------------

switch gcohparams.version
    case 'window'
        fprintf('Running gcoh on windows... ')
    case 'epoch'
        fprintf('Running gcoh on epoch(s)... ')
end

switch gcohparams.inference
    case 'none'
        disp('without inference.')
    case 'fdb'
        disp('with non-parametric FDB.')
end

%---------------------------------------
% BOOTSTRAP CROSS-SPECTRA 
% IF INFERENCE METHOD IS FDB, 
% ELSE MEAN ACROSS TAPERS
%---------------------------------------

if strcmp(gcohparams.inference,'fdb')
    
    fprintf('Computing FDB... ')
    
    bootPxy = nan([Nfreqs Nchans Nchans Nmw*Niter],'like',complex(0));
    
    for ff = 1:Nfreqs
        bootPxy(ff,:,:,:) = fdb_mtx(squeeze(Pxy(ff,:,:,:,:),Niter,gcohparams.version));
    end
    clear Pxy
    Pxy = bootPxy; clear bootPxy
    
    disp('Done.')
    
elseif strcmp(gcohparams.inference,'none')
    
    fprintf('Computing mean across tapers... ')
    
    Pxy = squeeze(mean(Pxy,2)); % becomes Nfreqs x Nchans x Nchans x Nmw
    
    disp('Done.')
    
end

%---------------------------------------
% RESHAPE CROSS-SPECTRA TO 3 DIMS
%---------------------------------------

if strcmp(gcohparams.inference,'none')
    
    Pxy_reshaped = zeros(Nchans,Nchans,Nfreqs*Nmw);
    c=1;
    for qq = 1:Nmw
        for loopfreq = frequencies
            Pxy_reshaped(:,:,c)=squeeze(Pxy(loopfreq,:,:,qq));
            c=c+1;
        end
    end
    clear Pxy;
    
elseif strcmp(gcohparams.inference,'fdb')
    
    Pxy_reshaped = zeros(Nchans,Nchans,Nfreqs*Nmw*Niter);
    c=1;
    for qq = 1:Nmw
        for loopfreq = frequencies
            Pxy_reshaped(:,:,c)=squeeze(Pxy(loopfreq,:,:,qq));
            c=c+1;
        end
    end
    clear Pxy;
    
end 

fprintf('Computing eigendecompositions... ')
%---------------------------------------
% RUN WINDOWED VERSION, NO FDB -- do first
%---------------------------------------

diags = zeros(Nchans, length(Pxy_reshaped));
eigvecs_r = zeros(Nchans,Nchans,length(Pxy_reshaped));
eigvecs_l = zeros(Nchans,Nchans,length(Pxy_reshaped));
nandat = find(squeeze(any(any(isnan(Pxy_reshaped),1),2))); % improve

for i = 1:size(Pxy_reshaped,3)
    
    if ismember(i,nandat)
        continue
    end
    
    [U,S,V] = eig(Pxy_reshaped(:,:,i)); % such that U*S*V' = Pxy
    [diags(:,i), svorder] = sort(diag(S),'descend');
    U = U(:,svorder);
    V = V(:,svorder);
    eigvecs_r(:,:,i) = U;
    eigvecs_l(:,:,i) = V;
    
end


eigenvalues = nan(Nmw,Nfreqs,Nchans);
eigenvectors_l = nan(Nmw,Nfreqs,Nchans,Nchans);
eigenvectors_r = nan(Nmw,Nfreqs,Nchans,Nchans);
c=1;
for qq = 1:Nmw
    for loopfreq = frequencies
        eigenvalues(qq,loopfreq,:) = diags(:,c);
        eigenvectors_l(qq,loopfreq,:,:)=eigvecs_l(:,:,c);
        eigenvectors_r(qq,loopfreq,:,:)=eigvecs_r(:,:,c);
        c=c+1;
    end
end

%---------------------------------------
% RUN WINDOWED VERSION WITH npFDB -- do third
%---------------------------------------



%---------------------------------------
% RUN EPOCHED VERSION, NO FDB -- do last or not at all
%---------------------------------------



%---------------------------------------
% RUN EPOCHED VERSION WITH npFDB -- do second
%---------------------------------------



disp('Done.')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gcohout.eigenvalues = eigenvalues;
gcohout.eigenvectors_r = eigenvectors_r;
gcohout.eigenvectors_l = eigenvectors_l;
gcohout.times = times;
gcohout.freqs = freqs;

%---------------------------------------
% END CODE
%---------------------------------------
end
