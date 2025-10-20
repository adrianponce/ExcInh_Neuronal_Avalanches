% Model's simulation 
%
% This code implements a Gillespie-based simulation of an extended stochastic 
% Wilson–Cowan model (Benayoun et al., 2010; de Candia et al., 2021). 
% Excitatory (E) and inhibitory (I) neurons interact through couplings wE and wI 
% with a nonlinear response function. The network architecture, 
% defined by experimental spatial coordinates and cell types, 
% used distance-dependent connectivity with exponentially decreasing probability. 
% Model spiking activity was converted into calcium-like transients 
% using a convolution model (Wei et al., 2020).
%
%
% This code uses:
%--------------------------
% Get_Connectivity_matrix.m
% Gillespie_EImodel.m
% SpikesToFluoresence.m
% Get_NonSpatialAvalanches.m
% plmle.m
% avalanche_analysis_EI.m
% avalanche_timing_EI.m
% neurons_type_and_coords.mat
%
% The function plmle.m belongs to MATLAB NNC (Neural Criticality and Complexity) Toolbox
% available here: http://www.nicholastimme.com/software.html
% (see also:  https://doi.org/10.3389/fphys.2016.00250)
%
%
% Adrian Ponce-Alvarez 11/07/2024
%--------------------------------------------------------------------------

% load example neurons' type and coordinates:
load('neurons_type_and_coords.mat','Type','xyz')

% Fixed model parameters:
%---------------------------------------------

% Inputs:
Io = 0.001;

% params:
beta_param = 1;
alpha_param = .1;
w0c = alpha_param/beta_param;

% transfer functions:
response_fn = @(x) beta_param*tanh(x).*(x>0);

% time window:
t_min = 0;
t_max = 1800;
Tran = 10; %transitory regime
n_batch = 1; %nb. of batches

% couplings:
wE = 7.1;
wI = 7; 
% length of connection probability 
lambda = 100;

% Construct connectivity:
%---------------------------------------------
[W,Ampli,NE,NI,typ,xyz] = Get_Connectivity_matrix(lambda,wE,wI,Type,xyz);
N = NE + NI;

% Convolution -> fluorescence signal
% (see Wei et al., 2020, PLoS Comput. Biol. 16(9): e1008198.):
%-------------------------------------------------------------
tauR = .5; % rise time (in the order of ms)
tauD = 3; % decay time (in the order of 0.5-2 s)
K = .6;  % between 0.5-3
q = 5;  % between 1-10
Fm = 10; % around 10
resol_F = 1/15;
            
% Run for reference point:
%-------------------------

I = Io*ones(N,1);    
    
% run simulation:
%---------------------------------------------
spike_times = [];
spike_ids   = [];
time = [];
Fds = [];
sh = 0;
for n = 1:n_batch
    % initial state:
       if n==1
       init_state = zeros(2,N);
       init_state(1,:) = rand(1,N)<.05;
       init_state(2,:) = ~init_state(1,:);
       else
       init_state = network_state;
       end
                
       fprintf('entering Gillespie algorithm... \n')
       [sp_times,sp_ids,network_state] = ...
                    Gillespie_EImodel(W,response_fn,beta_param,alpha_param,I,t_min,t_max,init_state);
            
       shift = (n-1)*t_max;
       spike_times = [spike_times (sp_times + shift)];
       spike_ids   = [spike_ids sp_ids];
                    
       % transform spikes to fluorescence signals:
       %-------------------------------------------------------------------
       [t,fds] = SpikesToFluoresence(sp_times,sp_ids,N,resol_F,tauR,tauD,K,q,Fm);
       if ~isempty(t)
          time = [time (t+sh)];
          sh = time(end);
          Fds = [Fds;fds];
          clear t fds
          fprintf('fluo done... \n')
       end       
end

%Thresholded calcium events:
threshold_fluo = 3*std(Fds);
%Binary events:
F = (Fds - repmat(mean(Fds),[size(Fds,1),1]));
% Binary events:
Raster_F = F > threshold_fluo;            
clear F Fds
            
% E/I ratio:
PopE = mean(Raster_F(:,1:NE),2);
PopI = mean(Raster_F(:,NE+1:end),2);
EIratio_F = nanmean( PopE./(PopE + PopI) );

% Avalanches:
X = sum(Raster_F,2);
th = floor(0.005*N);
[Size_F,Duration_F] = Get_NonSpatialAvalanches(X,th);
Nav_F = length(Size_F);


% <S>(T) function:
w = logspace(0,1.1*log10(max(Duration_F)),15);
Ts_F = nan(1,length(w)-1);
S_F = nan(1,length(w)-1);
for k=1:length(w)-1  
    ii=find(Duration_F>=w(k) & Duration_F<w(k+1));
    dur = Duration_F(ii);
    s   = Size_F(ii);
    Ts_F(k) = mean(dur);
    S_F(k) = mean(s);
end

 % Get exponents:
 %-----------------------------
 % Maximum Likelihood Estimation:
 tau_F = plmle(Size_F,'xmin',10,'xmax',10*N);
 Alpha_F = plmle(Duration_F,'xmin',min(Duration_F),'xmax',max(Duration_F));
 % Least-squares for <S>(T):
   cut = 1; %*resol; 
   X = Ts_F(Ts_F>=cut);
   Y = S_F(Ts_F>=cut);
   logx = log(X);
   logy = log(Y);
   if ~isempty(logx) && ~isempty(logy)
       p = polyfit(logx,logy,1);
       sigmaNuZ_F = 1/p(1);
   else
       sigmaNuZ_F = NaN;
   end

tau_exponent_F      = tau_F;
alpha_exponent_F    = Alpha_F;
sigmaNuZ_exponent_F = sigmaNuZ_F;


% Spatial definition (Ponce-Alvarez et al. 2018):
%----------------------------------------------------------------------
[Size_Space,Duration_Space,AvSize_E,AvSize_I,Avalanche,Av_cluster_sizes,Av_Center,Av_cluster_E,Av_cluster_I] = ...
    avalanche_analysis_EI(Raster_F,xyz,typ,10,3);

[Shapes,avg_shapes,Shapes_EI,avg_shapes_EI,AvD,numAv] = ...
    avalanche_timing_EI(Avalanche,Av_cluster_sizes,Av_cluster_E,Av_cluster_I);


if ~isempty(Size_Space)
% <S>(T) function:
w = logspace(0,1.1*log10(max(Duration_Space)),15);
Ts_Space = nan(1,length(w)-1);
S_Space = nan(1,length(w)-1);
for k=1:length(w)-1  
    ii=find(Duration_Space>=w(k) & Duration_Space<w(k+1));
    dur = Duration_Space(ii);
    s   = Size_Space(ii);
    Ts_Space(k) = mean(dur);
    S_Space(k) = mean(s);
end

% Get exponents:
%-----------------------------
% Maximum Likelihood Estimation:
tau_Space = plmle(Size_Space,'xmin',th,'xmax',10000);
Alpha_Space = plmle(Duration_Space*resol_F,'xmin',4*resol_F,'xmax',30);
% Least-squares for <S>(T):
cut = 4; %*resol; 
X = Ts_Space(Ts_Space>=cut);
Y = S_Space(Ts_Space>=cut);
logx = log(X);
logy = log(Y);
p = polyfit(logx,logy,1);
sigmaNuZ_Space = 1/p(1);
end
