% total number of neurons in network
N = 2;

% number of excitorty neurons in the network
Ne = 1;

% number of inhibitory neurons
Ni = N-Ne;

% random numbers for neurons
re = rand(Ne, 1);
ri = rand(Ni, 1);

% maximal delay in ms (time-window for STDP)
D = 20; 

% setting parameter values for neurons + randomization within some
% threshold
a = [0.02 * ones(Ne,1); 0.02+0.04*ri];
b = [0.2 * ones(Ne,1); 0.25-0.005*ri];
c = [-65+15*re.^2; -65*ones(Ni, 1)];
d = [8-6*re.^2; 2*ones(Ni,1)];

% the weight matrix holding placticity weights of connections
% between neurons
wmax=10;		% single synapse hard bound
%S = [0.5*rand(Ne+Ni,Ne), -rand(Ne+Ni,Ni)];
S = [wmax * ones(N, Ne), -wmax * ones(N,Ni)];



% delete neuron-neuron connections
for i=1:N
    S(i,i) = 0;
end;

% initial val  ues
v = -65 * ones(Ne+Ni, 1);
u = b.*v;

% we hold firing times using a two-column matrix. 
% First column describes spike-time and second which neuron fired.
tsteps = 1000;
t_grid = 1:tsteps;
v1_grid = zeros(1, tsteps);
v2_grid = zeros(1, tsteps);
dt_grid = linspace(-0.2,0.2,100);

for days = 1:100
    firings = [];
    
for t=t_grid
    % create random input external to the network
    I = [10; 0];
    
    % determine spiking neurons
    fired = find(v >= 30 );
    
    % add times of firings
    firings = [firings;t*ones(length(fired),1), fired];
    
    % reset fired neurons
    v(fired) = c(fired);
    u(fired) = u(fired) + d(fired);
    
    % add to the input I for each neuron a value = synaptic strengths of
    % all other neurons that fired in the last time
    I = I + sum(S(:,fired),2);
    
    % move simulation using Eulers method
    v = v + 0.5 * (0.04 * v.^2 + 5 .* v + 140 - u + I);
    v = v + 0.5 * (0.04 * v.^2 + 5 .* v + 140 - u + I);
    u = u + a .* (b .* v - u);
    
    v1_grid(t) = v(1);
    v2_grid(t) = v(2);
    
    % here we should change plasticity according to STDP rule
    S = 0.99999 * S;
    % differences in spiking times after stdp and neuron indices
    %stdp(firings(firings(:,2)==1 & (firings(:,1) >= t - D & firings(:,1) <= t + D ),1)-t);
    %firings(firings(:,2)==1 & (firings(:,1) >= t - D & firings(:,1) <= t + D ),2);
    %S = S + 
    %I = I + sum(S(:,fired),2);
end;

% normalizing
vplot(find(v1_grid >= 30)) = 30;
vplot(find(v2_grid >= 30)) = 30;

% plotting 
subplot(3,2,[1,2])
h1 = plot(t_grid,v1_grid);
%set(h,"linewidth", 2); 
axis([0 max(t_grid) -95 60])
xlabel('Time [ms]'); ylabel('Voltage [mv]');
title('Excitory neuron (1)');

subplot(3,2,[3,4])
h2 = plot(t_grid,v2_grid);
%set(h,"linewidth", 2); 
axis([0 max(t_grid) -95 60])
xlabel('Time [ms]'); ylabel('Voltage [mv]');
title('Inhibitory neuron (2)');

subplot(3,2,5); 
imagesc(S,[-wmax,wmax]); colormap(hot); colorbar
title('STDP-controlled weights'); xlabel('neuron index'); ylabel('neuron index');

subplot(3,2,6); 
cla('reset');
h3 = plot(dt_grid, stdp(dt_grid), '.');
title('STDP function'); xlabel('time difference dt [ms]'); ylabel('change in weight w_{ij}');

drawnow;
end;