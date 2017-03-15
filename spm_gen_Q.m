function [Q] = spm_gen_Q(P,X)
% Helper routine for spm_gen routines
% FORMAT [Q] = spm_gen_Q(P,X)
%
% P - parameters
% X - vector of between trial effects
% c - trial in question
%
% Q - trial or condition-specific parameters
%
% This routine computes the parameters of a DCM for a given trial, where
% trial-specific effects are deployed according to a design vector X. The
% parameterisation follows a standard naming protocol where, for example,
% X(1)*P.B{1} + X(2)*P.B{2}... adjusts P.A for all (input) effects encoded
% in P.B.
% P.BN and P.AN operate at NMDA receptors along extrinsic connections
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_Q.m 5922 2014-03-18 20:10:17Z karl $


% condition or trial specific parameters
%==========================================================================
Q  = P;


% % trial-specific effects on C (first effect only)
% %--------------------------------------------------------------------------
% try
%     Q.C = Q.C(:,:,1) + X(1)*P.C(:,:,2);
% end

% trial-specific effects on A (connections)
%--------------------------------------------------------------------------
for i = 1:length(X)
    
    % Alex: allow full trial specific inputs [e.g. roving]
    % trial-specific effects on C (first effect only)
    %--------------------------------------------------------------------------
    try
        if i > 1; Q.C(:,:,i) = Q.C(:,:,i) + X(i)*P.C(:,:,i-1); end
    end

    
    % extrinsic (driving) connections
    %----------------------------------------------------------------------
    for j = 1:length(Q.A)
        
        Q.A{j} = Q.A{j} + X(i)*P.B{i};
        
        % CMM-NMDA specific modulation on extrinsic NMDA connections
        %------------------------------------------------------------------
        try
            Q.AN{j} = Q.AN{j} + X(i)*P.BN{i};
        end
        
    end
    
    % modulatory connections
    %----------------------------------------------------------------------
    try
        Q.M  = Q.M + X(i)*P.N{i};
    end
    
    % intrinsic connections
    %----------------------------------------------------------------------
    try
        % Alex: vector specifying trial spec conns
        ACTI = [1 4 7 10]; ... SS-SS / II-II / SP-SP / DP-DP
       %Q.G(:,1)    = Q.G(:,1) + X(i)*diag(P.B{i});
        Q.G(:,ACTI) = Q.G(:,ACTI) +  repmat( X(i)*diag(P.B{i}) , [1 length(ACTI)] );
    end
    
end
