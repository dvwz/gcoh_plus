%MODALPROJ Modal projection
%Optional file header info (to give more details about the function than one line)
%Optional file header info (to give more details about the function than one line)
%
% Usage:  
%    [output1,output2] = modalproj(input1)
%    [output1,output2] = modalproj(input1,input2)
%
% Inputs:
%    mtspects - Multi-taper estimates (Ntimes x Nfreqs x Nchans)
%    eigvecs - Left eigenvectors (Nchans x Nchans x Nfreqs (x Niter))
%
% Outputs:
%    mdproj - Modal projection (0-1) (Ntimes x Nfreqs x Niter)
%    nwspect - Network-weighted spectrum (Ntimes x Nfreqs x Niter)
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Other requirements: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Copyright Apr-2020, David Zhou, dwzhou@mit.edu
% Last revision 06-Apr-2020
%------------------------------------------------
 
function [mdproj, nwspects] = modalproj(mtspects,eigvecs,varargin)

%---------------------------------------
% PROCESS INPUTS
%---------------------------------------
 
%VARARGIN
if nargin > 2
    modes=varargin{1};
else
    modes=[1];
end

Nfreqs = size(eigvecs,3);
Niter = size(eigvecs,4);
Ntimes = size(mtspects,1);
Nmodes = length(modes);
Nchans = size(mtspects,3);

% %VARARGOUT
% switch length(varargout)
%     case 1
%         execute_procedure1=1;
%     case 2
%         execute_procedure1=1;
%         execute_procedure2=1;
% end

%---------------------------------------
% MAIN BODY OF FUNCTION
%---------------------------------------
%% init output matrices

mdproj = nan(Ntimes,Nfreqs,Niter);
%mp_alpha = nan(Nfreqs,Ntimes,Niter);
trace_Pxys = nan(Nfreqs,Ntimes);
nwspects = nan(Ntimes,Nfreqs,Niter);

%% compute Pxy
J = mtspects; % Ntimes x Nfreqs x Nchans
Jx = conj(repmat(J,[1,1,1,Nchans]));
Jy = repmat(reshape(J,[Ntimes,Nfreqs,1,Nchans]),[1,1,Nchans,1]);
Pxys = Jx.*Jy; % Ntimes x Nfreqs x Nchans x Nchans
Pxys = permute(Pxys,[3 4 2 1]); % Nchans x Nchans x Nfreqs x Ntimes
mtspects = permute(mtspects,[1 3 2]); % Ntimes x Nchans x Nfreqs

%% pull eigenvectors

U_fs = eigvecs(:,modes,:,:); % Nchans x Nmodes x Nfreqs x Niter

for ti = 1:Ntimes
    
    % compute denominator from cross-spectrum trace
    trace_curr = zeros(Nfreqs,1);
    for chani = 1:Nchans
        trace_curr = trace_curr + squeeze(Pxys(chani,chani,:,ti));
        assert(length(trace_curr) == Nfreqs)
    end
    trace_curr = squeeze(trace_curr); % Nfreqs x 1
    assert(length(trace_curr) == Nfreqs)
    trace_Pxys(:,ti) = trace_curr;
    
    for ii = 1:Niter
%         % perform alpha-specific modal projection
% 
%         mp_curr = zeros(Nfreqs,1);
%         for fi = 1:Nfreqs
%             mp_curr(fi) = trace(U_alpha(:,:,ii)' * Pxys(:,:,fi,ti) * U_alpha(:,:,ii)); % (3 x Nchans) * (Nchans x Nchans) * (Nchans x 3) = (3 x 3)
%         end
%         mp_alpha(:,ti,ii) = real(mp_curr ./ trace_curr);

        % perform frequency-varying modal projection

        mp_curr = zeros(Nfreqs,1);
        for fi = 1:Nfreqs
            mp_curr(fi) = trace(U_fs(:,:,fi,ii)' * Pxys(:,:,fi,ti) * U_fs(:,:,fi,ii)); % (Nmodes x Nchans) * (Nchans x Nchans) * (Nchans x Nmodes) = (Nmodes x Nmodes)
        end
        mdproj(ti,:,ii) = real(mp_curr ./ trace_curr);

        % perform network-weighted spectrograms
        for fi = 1:Nfreqs
            nweight = diag(U_fs(:,:,fi,ii) * U_fs(:,:,fi,ii)'); % (Nchans x Nmodes) x (Nmodes x Nchans)
            nwspects(ti,fi,ii) = real(diag(Pxys(:,:,fi,ti))' * nweight);
        end

    end
end

% % median over iterations
% mdproj = nanmedian(mdproj,3)';
% %mp_alpha = nanmedian(mp_alpha,3)';
% nwspects = nanmedian(nwspects,3)';

%---------------------------------------
% SUBFUNCTIONS
%---------------------------------------


%---------------------------------------
% END CODE
%---------------------------------------
end
