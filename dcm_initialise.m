function D = dcm_initialise()
% wrapper for ERP_DCM_AS for batching big jobs with several groups and
% models
%
% must set up MODELS.m        % model spaces for n models
%         and GroupDataLocs.m % data locations & info for n groups
%
% AS



% check MODELS.m and GroupDataLocs.m are set!
H = pwd; 
M = which('MODELS.m');      if isempty(M); disp('Check MODELS.m in path'); return; end
G = which('GroupDataLocs'); if isempty(G); disp('Check GroupDataLocs.m in path'); return; end


% get what we need
M      = MODELS;        m = length(M);
[gf,p] = GroupDataLocs; g = length(p);



% ERP_DCM_AS(m, StNo,G) where: 
% m    = the models you want to run (int or vector)
% StNo = subjects
% G    =  Group (corresponding to GroupDataLocs.m)

for i = 1:length(gf) % groups
    
    nsub = length(dir([gf{i} '/' p(i).f]));
    
    ERP_DCM_AS(1:m,1:nsub,i);
    
end