function [spike_times,spike_ids,network_state] = Gillespie_EImodel(W,response_fn,beta,alpha,I,t_min,t_max,init_state)

% This function is based on Edward Wallace's code, publicly available here:
% https://github.com/ewallace/stochsimcode
%
% This code simulates a stochastic Wilson-Cowan model with the Gillespie algorithm.
% The connectivity of E and I neurons is given by matrix W (N-by-N, where N = NE + NI). 
% The background inputs are given by I (N-dimensional). 
% State transition rates are given by alpha and the tranfert function (response_fn, handle function)
% 
%
% Inputs:
% -W: connectivity matrix, N-by-N: W(i,j) is synaptic
% weight from the jth neuron to the ith neuron.
% -response_fn: handle for response function (sigmoid, hyperbolic tan, etc.)
% -I: backgound inputs, N-dimensional array
% -alpha:  rate at which active neurons decay to quiescent state, N-dimensional array
% -beta: is the height of the response function, N-dimensional array
% -init_state: is the initial state vector, 2-by-N
%
% Outputs:
% - spike_times: spike times
% - spike_ids: neuron IDs associated to spike times
% - network_state: last network state (useful for batches)
%
% References:
% Citation: Benayoun M, Cowan JD, van Drongelen W, Wallace E (2010) 
% Avalanches in a Stochastic Model of Spiking Neurons. 
% PLoS Comput Biol 6(7): e1000846. 
% https://doi.org/10.1371/journal.pcbi.1000846
%
% Wallace E, Benayoun M, van Drongelen W, Cowan JD. 
% Emergent Oscillations in Networks of Stochastic Spiking Neurons 
% PLoS ONE, May 6, 2011. DOI: 10.1371/journal.pone.0014804
%--------------------------------------------------------------------------

N = size(W,1);

% Calculate expected number of events
factor=10;
expected_events=N*(t_max-t_min)*factor;

% Initialize vectors for update times, label of updated neuron,
% and new state of updated neuron
times = zeros(1, expected_events);
updates = times;
new_states = zeros(2, expected_events);

% Set event counter to 0 and simulation time to initial time.
event_no = 0;
curr_time = t_min;
dt = 0;

% initialize network state vector - we'll keep one vector
% for the active neurons and another for the quiescent ones
network_state = init_state;
active = init_state(1,:)';
quiescent = init_state(2,:)';

% Calculate vector of transition rates at initial time
currents = W*active + I;
trans = beta .* (active==0) .*feval(response_fn,currents) + ...
    alpha.*(active==1);
cum_trans=cumsum(trans);


% Main loop: update according to Gillespie algorithm, with rates
% specified by trans, until time t_max is exceeded.

while (curr_time < t_max)
    curr_time = curr_time + dt;
    
    % Call gillespie to pick update time, neuron updated, and new state

    % Calculates total network transition rate, as sum of transition
    % rates of all neurons, i.e. last element of cum_trans.
    total_trans = cum_trans(end);

    % timestep is exponential R.V. with parameter total_trans
    dt = -log(rand)/total_trans;

    % pick random variable uniform on (0, total_trans), and select
    %  neuron with number i_update as least i with
    % test_variable < cum_trans(i)
    test_variable = total_trans*rand;

    i_update = find(cum_trans >= test_variable,1,'first');

    % Pick new state for neuron i_update
    new_state = double(~network_state(:,i_update));

    active(i_update) = new_state(1);
    quiescent(i_update) = new_state(2);
    
    % change transition rates of neurons affected by spike
    if(new_state(1) == 1)
        currents = currents + W(:,i_update);
    elseif(new_state(1)==0)
        currents = currents - W(:,i_update);
    end
    trans = beta .* (active==0) .*feval(response_fn,currents) + ...
        alpha.*(active==1);
    
    
    % calculate cumulative sum of transition probabities
    cum_trans=cumsum(trans);
    
    
    event_no=event_no+1;
    times(event_no)=curr_time+dt;
    updates(event_no) = i_update;
    new_states(:, event_no)=new_state;
    network_state = [ active' ; quiescent' ];
    
    % If network can make no further transitions, end simulation.
    if(cum_trans(N)==0)
        break
    end
end
% End of main simulation loop

event_no = event_no - 1;
times=times(1:event_no);
updates=updates(1:event_no);
new_states=new_states(:, 1:event_no);


% Takes vectors of all times, updates, and new states, outputs 
% reduced vector containing only spike times and info, i.e. 
% transitions from quiescent to active.
spikes = find(new_states(1,:));

spike_times = times(spikes);
spike_ids = updates(spikes);


return
