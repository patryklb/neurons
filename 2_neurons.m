% Model of the network consisting of two Leaky-integrate-and-fire neurons(LIF) 
% with noise and refractory current.

N = 2;                                          % number of neurons
T = 6;                                        % time of simulation
dt = 0.01;                                      % time step
tgrid = 0:dt:T+dt;                              % time grid
n = numel(tgrid);                               % number of time steps
t_rest = 0;                                     % initial refractory time
spikes = zeros(N,n);                            % stores data about spikes

% LIF parameters:
Vm = zeros(N,n);                                % potential over time
Cm = [5;5];                                  	% capacitance in [uF]
Rm = [5;5];                                     % resistance in [kOhm]
tau_m = Rm.*Cm;                                 % time constant (msec)
tau_ref = 4;                                    % refractory period
Vth = [1.6;1.6];                                    % spike threshold
I = [10*ones(1,n); 10*ones(1,n)];                % input current (A)
V_spike = 0.5;    
V_cut = [0;0];

Vm(:,1) = 1;

% TIME DYNAMICS LOOP

for i=1:n-1
  Vm(:,i+1) = Vm(:,i) + (-Vm(:,i)+I(:,i).*Rm./tau_m)*dt; 
  
  if(Vm(1,i+1)>Vth(1,1))
    Vm(1,i+1) = V_cut(1,1);  
    I(2, i+1) =  I(2, i+1) + Vth(1,1)/Rm(2);
  end
  if(Vm(2,i+1)>Vth(2,1))
    Vm(2,i+1) = V_cut(2,1);
    I(1, i+1) =  I(1, i+1) + Vth(2,1)/Rm(2);
  end
end
subplot(2,1,1)
h1 = plot(tgrid,Vm(1,:));
set(h1,"linewidth", 2); 
xlabel('Time [ms]'); ylabel('Voltage [mv]');
title('Neuron 1');
subplot(2,1,2)
h2 = plot(tgrid,Vm(2,:));
set(h2,"linewidth", 2);
xlabel('Time [ms]'); ylabel('Voltage [mv]');
title('Neuron 2');