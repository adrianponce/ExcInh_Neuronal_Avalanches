function [W,Ampli,NE,NI,typ,xyz] = Get_Connectivity_matrix(lambda,wE,wI,Type,xyz)

% Constructs connectivity
% Inputs:
% -lambda param.: length of connection prob. decay
% -wE and WI: exc. and inh. couplings
% -Type: neurons' type (-1 = Ch; 0 = Inh; 1 = Exc)
% -xyz: neurons' coordinates
%
% Outputs:
% - W : connectivity matrix
% - Ampli : balanced amplification index
% - NE, NI : nb. of E and I neurons
% - type: cell-type array
%
%
% Adrian Ponce-Alvarez 11/07/2024
%--------------------------------------------------------------------------

% Spatially embeded connectivity

    
    XYZ = xyz(Type>-1,:);
    typ = Type(Type>-1);
    
    N = length(typ);
    NE = sum(typ==1);
    NI = sum(typ==0);
        
        D = zeros(N);
        for i =1:N
            ri = XYZ(i,:);
            for j=1:N
                rj = XYZ(j,:);
                D(i,j) = sqrt( sum((ri-rj).^2) ); % Euclidean dist.
            end
        end


    WEE = rand(NE) < exp(-D(1:NE,1:NE)/lambda);
    WIE = rand(NI,NE) < exp(-D(NE+1:end,1:NE)/lambda);
    WEI = rand(NE,NI) < exp(-D(1:NE,NE+1:end)/lambda);
    WII = rand(NI) < exp(-D(NE+1:end,NE+1:end)/lambda);

    WEE(1:NE+1:end) = 0; WII(1:NI+1:end) = 0;
    
    
    WEE_norm = WEE./repmat(sum(WEE,2)/wE,[1,size(WEE,2)]);
    WEI_norm = WEI./repmat(sum(WEI,2)/wI,[1,size(WEI,2)]);
    WIE_norm = WIE./repmat(sum(WIE,2)/wE,[1,size(WIE,2)]);
    WII_norm = WII./repmat(sum(WII,2)/wI,[1,size(WII,2)]);
    
    W = [WEE_norm -WEI_norm;WIE_norm -WII_norm];
    W = single(W);
    
    T = schur(W);
    [~,d] = eig(T);
    d = diag(d);
    Y = T*T';
    Ampli = 1 - sum( abs(d).^2 )/trace(Y);
    xyz = XYZ;


return