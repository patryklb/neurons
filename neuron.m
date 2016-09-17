% The Leaky-integrate-and-fire neuron model (LIF) 
% with noise with refractory current.
 
function [Vm,tgrid,spikes] = neuron(I,Rm,Cm,Vth)
%Rm = 1;			% resistance in [kOhm]
%Cm = 1;			% capacitance in [uF]
%Vth = 1;			% spike threshold

T = 30;				% 
dt = 0.01;			% time step
tgrid = 0:dt:T;			% time grid 
t_rest = 0;			% initial refractory time

% LIF parameters:
Vm = zeros(1,length(tgrid));	% potential over time

tau_m = Rm*Cm;			% time constant (msec)
tau_ref = 5;			% refractory period (msec)
V_spike = 0.5;			% spike delta
N = numel(tgrid);		% number of time-steps
spikes = 0;

% MAIN TIME LOOP
for i=1:N-1
  if tgrid(i)>t_rest
    Vm(i+1) = Vm(i) + (-Vm(i)+I*Rm) / tau_m * dt;
  end
  if (Vm(i+1) >= Vth)
    Vm(i+1) += V_spike;
    t_rest = tgrid(i) + tau_ref;
    spikes += 1;
  end
end


% plot membrane potential Vm
h = plot(tgrid, Vm) 1
set(h,"linewidth", 2); 
title('Leaky Integrate-and-Fire model');
ylabel('Voltage [V]'); xlabel('Time (msec)');
ylim([0,2]);
end

