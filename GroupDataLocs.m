function [s,p] = GroupDataLocs
% For use with ERP_DCM_ALEX.m
% This func for finding your data (over groups).
%
% AS2016

s{1} = '/Users/Alex/Desktop/Roving/NewERP'; % group folder
p(1).xU.X     = eye(8);            % roving
p(1).xU.name  = {'Deviant' 'r1','2','3','4','5','6','7'}; % condition names
p(1).tCode    = 1:8;             % condition codes in SPM
p(1).d        = '*';            % subj directories
p(1).f        = 'ALL_Six*.mat';    % file to find



% s{1} = '/imaging/as08/FTD_MM_DCM_moresubs/cons1/'; % group folder
% p(1).xU.X     = [1 0; 0 1];        % deviant on-off
% p(1).xU.name  = {'Deviant' 'Std'}; % condition names
% p(1).tCode    = [1 2];             % condition codes in SPM
% p(1).d        = 'meg*';            % subj directories
% p(1).f        = 'All_Six*.mat';    % file to find
% 
% 
% s{2} = '/imaging/as08/FTD_MM_DCM_moresubs/cons2/';
% p(2).xU.X     = [1 0; 0 1]; % design matrix
% p(2).xU.name  = {'Deviant','Std'};
% p(2).tCode    = [5 4]; % for cons2: 4 = FreqStd and 5 = 'Freq'
% p(2).d        = 'meg*';            % subj directories
% p(2).f        = 'All_Six*.mat';
% 
% s{3} = '/imaging/as08/FTD_MM_DCM_moresubs/cons3/';
% p(3).xU.X    = [1 0; 0 1];
% p(3).xU.name = {'Deviant','Std'};
% p(3).tCode   = [11 10]; % 11=Freq,10=FreqStd
% p(3).d        = 'meg*';            % subj directories
% p(3).f        = 'All_Six*.mat';
% 
% 
% 
% 
% s{4} = '/imaging/as08/FTD_MM_DCM_moresubs/ftd1/';
% p(4).xU.X    = [1 0; 0 1];
% p(4).xU.name = {'Deviant','Std'};
% p(4).tCode   = [1 2]; % simple mismatch: deviant vs std
% p(4).d        = 'meg*';            % subj directories
% p(4).f        = 'All_Six*.mat';
% 
% s{5} = '/imaging/as08/FTD_MM_DCM_moresubs/ftd2/';
% p(5).xU.X    = [1 0; 0 1];
% p(5).xU.name = {'Deviant','Std'};
% p(5).tCode   = [1 2];
% p(5).d        = 'meg*';            % subj directories
% p(5).f        = 'All_Six*.mat';
% 
% s{6} = '/imaging/as08/FTD_MM_DCM_moresubs/ftd3/';
% p(6).xU.X    = [1 0; 0 1];
% p(6).xU.name = {'Deviant','Std'};
% p(6).tCode   = [9 8];
% p(6).d        = 'meg*';            % subj directories
% p(6).f        = 'All_Six*.mat';
% 
% s{7} = '/imaging/as08/FTD_MM_DCM_moresubs/ftd4/';
% p(7).xU.X    = [1 0; 0 1];
% p(7).xU.name = {'Deviant','Std'};
% p(7).tCode   = [2 7];
% p(7).d        = 'meg*';            % subj directories
% p(7).f        = 'All_Six*.mat';
% 
% % data from pilot
% s{8} = '/imaging/as08/FTD_MM_DCM_pilot/FTD/';
% p(8).xU.X    = [1 0; 0 1];
% p(8).xU.name = {'Deviant','Std'};
% p(8).tCode   = [2 1];
% p(8).d        = 'meg*';            % subj directories
% p(8).f        = 'All_Six*.mat';
