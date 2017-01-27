function ERP_DCM_AS(m,StNo,G)
% Top level function for ERP DCM. Calls a  number of updated SPM function,
% stored in /imaging/as08/Roving/NEW_DCM_ERP_SPM12/cons1/. These include
% hacked version of spm_dcm_erp, spm_dcm_erp_data & spm_lx_erp, so add them
% to your paths.
%
% This script is accompannied by 2 further functions:
% 1) MODELS.m, which contains all the model architectures you want to run
% 2) GroupDataLocs.m, which contains the locations of the data you want to
% run. This can handle multi-group data.
%
% Also note that you will want to scroll down and look at:
% 1) CustomPriors - here are some custom priors [currently for CMC model]
% 2) PrepData     - IMPORTANT: options for the empirical data
% 3) checktrialcodes - special routine for checking trial codes
%
% USAGE:
% ERP_DCM_AS(m, StNo,G) where: 
% m    = the models you want to run (int or vector)
% StNo = subjects
% G    =  Group (corresponding to GroupDataLocs.m)
% 
% e.g. ERP_DCM_AS(1:21,1:16,1); will run models 1:21, subs 1:16 in
% group1.
%
% AS2016 [DCM]


% ensure SPM12 + custom functions
try ls /imaging ; CBU = 1; clc; catch CBU = 0; end


if CBU; addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest'));
        addpath('/imaging/as08/Roving/NEW_DCM_ERP_SPM12/');
        addpath('/imaging/as08/Roving/NEW_DCM_ERP_SPM12/cons1/');
else    addpath(genpath('~/Downloads/spm12'));
        %addpath('~/Desktop/Roving/NewERP/');
        addpath('/Volumes/Extra/Roving/NewERP/');
end


Mst       = m(end);         % Stop at model m(end)
[M,C,L,B] = MODELS;         % Model spaces
[s,p]     = GroupDataLocs;  % Data locations

s = s{G};
p = p(G);


% IGNORE... [calls subs functions]
lM         = length(M);
try if Mst ~=0;
       lM  = Mst;
    end;
end
for i = m(1):lM
    
    fprintf('\nrunning model %d\n',i);
    try kount = kount + 1; catch kount = 1; end
    if  kount > 1; StNo(1) = 1; end
    Do(M{i},C{i},L{i},B{i},i,s,p,StNo,p.d,p.f);
    
end

end

function DCM = PrepData(DCM,Ns,tCode)
% Sets options for the empirical data, with new options including
% baselining, filtering and using a set % of trials.

DCM.M.U            = sparse(diag(ones(Ns,1)));  ... ignore [modes]
DCM.options.trials = tCode;                     ... trial code [GroupDataLocs]
DCM.options.Tdcm   = [1 350];       ... peristimulus time
DCM.options.Fdcm   = [4 48];        ... frequency window
DCM.options.D      = 2;             ... downsample
DCM.options.han    = 1;             ... apply hanning window
DCM.options.h      = 4;             ... can't remember
DCM.options.DoData = 1;             ... leave on [custom]
DCM.options.Bdcm   = [-200 0];      ... baseline times [new!]
DCM.options.Fltdcm = [1 15];        ... bp filter [new!]

DCM.options.analysis      = 'ERP';  ... analyse type
DCM.xY.modality           = 'LFP';  ... ECD or LFP data? [LFP]
DCM.options.spatial       = 'LFP';  ... spatial model [LFP]
DCM.options.model         = 'CMC';  ... neural model
DCM.options.Nmodes        = length(DCM.M.U);

DCM.xY.name     = {'L V1','R V1','L pSTS', 'R FFA'};
DCM.Sname       = DCM.xY.name';
DCM.xY.Ic       = [1:Ns];

ForceCheck = 0; % this will force checktrialcodes [see below before using]

if ForceCheck; 
    DCM = checktrialcodes(DCM);
end

DCM = spm_dcm_erp_data(DCM);


end

function Do(M,C,L,B,nm,s,p,StNo,d,f)
% Make the DCM structure. 
% This will contain empirical data in DCM.xY, experiment [design] info in 
% DCM.xU, function handles to the generative [DCM.M.f], forward [DCM.M.G/g], 
% integrator [DCM.M.IS] and feature selection [DCM.M.FS] as well as the
% prior estimate and variances for those functions in DCM.M.pE/gE/pC/gC.


% IGNORE ME
h     = pwd; s = s; cd(s); % [for save]
DDir  = dir(d); DDir = {DDir.name};
for n = 1:length(DDir); try fi{n} = [s ls(deblank([DDir{n} '/' f]))]; end;end; N = n; cd(h);
Data  = fi;


% Start from subject n, or 1:
if isempty(StNo) || exist('StNo') == 0
     ns = 1;
else ns = StNo;
     fprintf('starting from subject %d',ns(1));
end


    % loop over subjects
    %----------------------------------------------------------------------
    for s = ns(1):ns(end)
    fprintf('\nrunning subject %d of %d\n',s,length(Data));
    
    
        % data naming & desing matrix
        %--------------------------------------------
        DCM          = [];
        [fp fn fe]   = fileparts(Data{s});
        DCM.name     = genvarname(['testERP_Mod_',num2str(nm),'_', fn, fe]);
        DCM.name     = [fp '/' DCM.name];
        
        if exist([DCM.name '.mat']) == 2; 
            fprintf('found = skipping sub %d',s);
        else
            
        DCM.xY.Dfile = Data{s};
        Ns           = 4;
        DCM.xU.X     = p.xU.X;       ... design matrix
        DCM.xU.name  = p.xU.name;    ... condition names
        tCode        = p.tCode;      ... condition index (in SPM)
        
        
        % Extrinsic connectivity - model spaces
        %---------------------------------------------
        ALL      = M;                          ... full connectivity
        DCM.B    = {B};                        ... trial specific
        DCM.A{3} = L;                          ... lateral [modulatory]
        DCM.C    = C;                          ... [exogenous] inputs
        
        %DCM.A{1} = triu(ALL);                  ... forward
        %DCM.A{2} = tril(ALL);                  ... backward
        DCM.A{1} = (ALL==1);
        DCM.A{2} = (ALL==2);
        
        DCM.B(2:length(DCM.xU.X))=DCM.B;
        
        % Functions
        DCM = CustomPriors(DCM,Ns);    % anything non-default / built in
        DCM = PrepData(DCM,Ns,tCode);  % evaluate empirical data
        DCM = MakeAndInvert(DCM);      % invert the model & save it
        end
    end

end

function DCM = CustomPriors(DCM,Ns)

% This will put a bunch of custom (non-DCM/SPM) functions or priors into a
% new structure in the DCM. These will be copied over in the call to
% spm_dcm_csd. As such this will require a custom version of spm_dcm_csd
% which includes a line to detect whether this struct has been provided and
% to copy over it's contents prior to inversion.

DCM.CUSTOM      = [];
DCM.CUSTOM.f    = 'spm_fx_cmcKRISH13';          ... generative model  (.m)
DCM.CUSTOM.pE   = [];                           ... intrinsic priors
DCM.CUSTOM.pE.G = zeros(Ns,13);                 ... [local coupling]
DCM.CUSTOM.pC.G = zeros(Ns,13);                 ... variance [off]

Self = find(diag(speye(Ns).*( DCM.A{1}+DCM.A{2} ))); % SP gain only if in model
DCM.CUSTOM.pC.G(Self,7) = 1/8;

% Remove self conns from extrinsic connectivity matrix?
DCM.A{1} = DCM.A{1}.*~eye(Ns);
DCM.A{2} = DCM.A{2}.*~eye(Ns);

DCM.CUSTOM.pE.T = zeros(Ns,4);                  ... population time const
DCM.CUSTOM.pC.T = zeros(Ns,4)+1/8;              ... variances
DCM.CUSTOM.gE.J = sparse([1 3 7],1,[.2 .8 .2],8,1)'; ... contributing states
DCM.CUSTOM.gC.J = sparse([1 3 7],1,1/8       ,8,1)'; ... variances
%DCM.CUSTOM.gC.J = repmat(DCM.CUSTOM.gC.J,[3 1]);     ... contrib sources
%DCM.CUSTOM.gE.J = repmat(DCM.CUSTOM.gE.J,[3 1]);     ... variance


end



function DCM = checktrialcodes(DCM)
% This will load up an SPM file and grab the condition names. It will check
% that your suplied codes exist and match the names. Provide strings in 'd' and 'S'
% below (which were for standard and deviant in my mismatch data example),
% to auto-select the correct codes if yours are wrong.

D = spm_eeg_load(DCM.xY.Dfile);
C = D.condlist;
n = DCM.options.trials;

%d = {'Deviant','Dev','FreqDev','FreDev','Freq'};  % possible names
%S = {'Standard','FreqStd','FreStd','StdF','Std'}; % possible names
d = {'Deviant','Dev','DevAll'};
S = {'Standard','Std','StdAll'};

try C(n(1)); catch [n(1),DCM.options.trials(1)] = deal(1); end
try C(n(2)); catch [n(2),DCM.options.trials(2)] = deal(1); end

% SEARCH FOR MATCHING TRIAL NAME: DEVIANT
if ~ismember(C(n(1)),d);
    for i = 1:length(d)
        f = (strfind(C,d{i}));
        F = find(~cellfun(@isempty,f));
        if ~isempty(F); 
            DCM.options.trials(1) = F; 
            fprintf('\nchanging trial code for deviant to %d\n',F);
            DoGrabAll = 0;
            break;
        end
    end
    if isempty(F); 
        DoGrabAll = 1; 
    end
else DoGrabAll = 0;
end

% IF DEVIANT NOT FOUND, FIND DEVIANTS SEPERATELY
if DoGrabAll == 1
    All = {'Dur','Freq','Gap','Int','Side'};
    if sum(ismember(C,All)) > 4
        %DCM.options.trials(1)
        DCM.options.TRIALZ1 = find(ismember(C,All));
    end
end

clear DoGrabAll
clear F f

% SEARCH FOR MATCHING TRIAL NAME: DEVIANT
if ~ismember(C(n(2)),S);
    for i = 1:length(S)
        f = (strfind(C,S{i}));
        F = find(~cellfun(@isempty,f));
        if ~isempty(F); 
            if length(F) > 1; F = find(strcmp(C,'StdAll')); end
            DCM.options.trials(2) = F; 
            fprintf('\nchanging trial code for standard to %d\n',F);
            DoGrabAll = 0;
            break; 
        end
    end
    if isempty(F)
        DoGrabAll = 1;
    end
else DoGrabAll = 0;
end

% IF STD NOT FOUND, FIND STDS SEPERATELY
if DoGrabAll == 1
    All = {'StdD','StdF','StdG','StdI','StdS'};
    if sum(ismember(C,All)) > 4
        %DCM.options.trials(2) 
        DCM.options.TRIALZ2 = find(ismember(C,All));
    end
end
   
end

function DCM = MakeAndInvert(DCM)

DCM.CUSTOM.nograph   = 1;
DCM = spm_dcm_erp(DCM);

end