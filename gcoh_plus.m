%GCOH_PLUS Top-level function of the gcoh2 toolbox
% Description:
%    Could be a wrapper function, or object-oriented class method
%    Desired functionality: creates output data path, calls all versions of
%    the toolbox (per-window/per-epoch, single/bootstrap/null mtxspctrms,
%    parallel/nonparallel,...)
%    Settings: # of components to compute, multitaper settings
%    
% Usage:  
%    [output1,output2] = gcoh_plus(input1)
%    [output1,output2] = gcoh_plus(input1,input2)
%
% Inputs:
%    data - time series (in form nsamples x nchannels)
%    movingwin - (in the form [window winstep])
%    epoch - (in the form [tstart tend], or a nsegs x 2 array, or [])
%    params - structure with fields tapers, Fs, fpass, ...
%    gcohparams (optional) - structure with fields savedir, outputdir,
%        infmethod ('none', 'bootstrap', 'null'), parallel (boolean),
%        niter, 
%
% Outputs:
%    eigvecs - left eigenvectors (in form nchannels x nchannels)
%    eigvals - eigenvalues (in form nchannels x 1)
%    gcohmat - global coherence spectrum or time-frequency plot
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
 
function [outputs,varargout] = gcoh_plus(data, ...
                                         movingwin, ...
                                         epoch, ...
                                         params, ...
                                         gcohparams)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROCESS INPUTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% PARAMS
if ~isstruct(params)
    error('gcoh2:paramsstruct','Error:  params must be a structure.')
elseif ~isfield(params, 'Fs')
    error('gcoh2:missingFs','Error:  Please enter a value for params.Fs')
end

% default multitaper parameters for Chronux
if ~isfield(params, 'pad'); params.pad = 0; end
if ~isfield(params, 'tapers'); params.tapers = [3 2]; end
if ~isfield(params, 'fpass'); params.fpass = [0,(Fs/2)]; end
if ~isfield(params, 'err'); params.err = 0; end
if ~isfield(params, 'trialave'); params.trialave = 0; end

% GCOHPARAMS
if ~isstruct(gcohparams)
    error('gcoh2:gcohparamsstruct','Error:  gcohparams must be a structure.')
end

if ~isfield(gcohparams, 'savedir'); gcohparams.savedir = pwd; end
if ~isfield(gcohparams, 'outputdir'); gcohparams.outputdir = pwd; end
if ~isfield(gcohparams, 'infmethod'); gcohparams.infmethod = 'none'; end
if ~isfield(gcohparams, 'niter')
    switch gcohparams.infmethod
        case 'none'
            gcohparams.niter = 1;
        case 'bootstrap'
            gcohparams.niter = 200;
        case 'null'
            gcohparams.niter = 200;
    end
end
if ~isfield(gcohparams, 'parallel'); gcohparams.parallel = 0; end

% %INPUT ERRORS
% if inputerror1
%     error('gcoh2:inputError','Error:  Message goes here.');
% end
% 
% %GLOBAL VARIABLES
% globalvar1=some_value;
% globalvar2=some_other_value;
% 
% %OTHER INPUT PROCESSING
%     %insert here


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try 
    %PROCEDURE 1 NAME
    %code goes here
catch err
        % issue a warning
    warning('gcoh2:warningType',...
        ['Warning: The following went wrong in Procedure 1 Name --\n' err.message ]);
    %error handling code here

        % or issue an error & exit function:
    error('gcoh2:errorType',...
        ['Error: The following went wrong in Procedure 1 Name --\n' err.message ]);
end

try
    %PROCEDURE 2 NAME
    %code goes here
catch err
        % issue a warning
    warning('gcoh2:warningType',...
        ['Warning: The following went wrong in Procedure 2 Name --\n' err.message ]);
    %error handling code here

        % or issue an error & exit function:
    error('gcoh2:errorType',...
        ['Error: The following went wrong in Procedure 2 Name --\n' err.message ]);
end


%---------------------------------------
% SUBFUNCTIONS
%---------------------------------------


%---------------------------------------
% END CODE
%---------------------------------------
end
