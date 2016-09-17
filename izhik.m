function izhik(a, b, c, d)
%%%%%%%%%%%%%%% Izhikevich neuron model %%%%%%%%%%%%%%%%%%%%%% 
figure(1)
% For RS neuron set: a=0.02; b=0.2;  c=-65;  d=6;
v = -70;    u = b*v;
vv = [];    uu = [];    ii = [];
tau = 0.25;
tgrid = 0:tau:100;
T = tgrid(end)/10;

% defining current injection
ii = [zeros(1,T) 14*ones(1,tgrid(end)-T)];

% starting time loop (dynamics)
for t=tgrid
    %%%%%%%%%%%%%%%%%%%%%%%%%CURRENT INJECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = ii(end);
    % dynamical equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = v + 0.04*v*v+5*v+140-u+I;
    u = u + a*(b*v-u);
    % spike condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (v > 30)
        vv(end+1) = 30;
        v = c;
        u = u+d;
    else
        vv(end+1) = v;
    end
    uu(end+1) = u;
end;
h = plot(tgrid,vv, 1:100, -90+ii, 'LineWidth',2)
%set(h,"linewidth", 2); 
axis([0 max(tgrid) -95 60])
 xlabel('Time [ms]'); ylabel('Voltage [mv]');
 title('Itzhikevitz neuron model');

end
 