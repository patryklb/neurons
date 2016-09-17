% network of 2 IF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETETRS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1 = 0; deli = 1.1; trefrac = 0.25; dtr = 0.1; taus = 1.6; tauf = 0.1;
gsyn1 = 0.08; gsyn2 = 0.08; vsyn = 0.5;
gel = 0; bspike = 0.1;
vreset = 0; iall = 0.92; vtresh = 1;
iapp1 = 0.1; iapp2 = 0.1;
istim = 1; pon = 10; deltat = 0.03; inoise = 0; ponoise = 0; dnoise = 100;
rnoise = 0; signoise = 0.2; d = 0.1; iscale = 1; dd = 0.1;
save('params.mat')

%% dd - delay parameter
%% dt - real time-step
dd = 1;
dt = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINING ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tsteps = 1000;
t_grid = 1:tsteps;

v = zeros(2,tsteps);
tr = zeros(2,1);
s = zeros(2,tsteps);
tt = 0;

v(1,1) = 0;
v(2,1) = 0.2;
tr(1) = -5;
tr(2) = -5;
s(1,1) = 0.011;
s(2,1) = 0.011;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remember - delayed s1, s0
% time is not ok!
for t = t_grid-1
    if t < dd
        t = dd;
    end;
    
    % spike conditions
    if (v(1, t) - vtresh) >= 0  
        disp('1: trach!')
        v(1, t) = vreset;
        %tr(1) = tt;
    end;
    
    if (v(2, t) - vtresh) >= 0  
        disp('2: trach!')
        v(2, t) = vreset;
        %tr(2) = tt;
    end;
        
    
    % Euler step (double step for numerical stability)
    v(1,t+1) = v(1,t) + 0.5 * dt * ((-v(1,t) - gsyn1 * s(2,t-dd+1) * (v(1,t) - vsyn) + iall + iapp1) * heav(tt-(tr(1)+trefrac)) + gel * ((v(2,t) - v(1,t)) + bspike * heav((tr(2)+dtr)-tt)) + normrnd(rnoise, signoise) * iappp(dt*t,iscale*rand) + iapp(dt*t, istim));
    v(2,t+1) = v(2,t) + 0.5 * dt * ((-v(2,t) - gsyn2 * s(1,t-dd+1) * (v(2,t) - vsyn) + iall + iapp2) * heav(tt-(tr(2)+trefrac)) + gel * ((v(1,t) - v(2,t)) + bspike * heav((tr(1) + dtr)-tt)));
    % + normrnd(rnoise, signoise) * iappp(dt*t, iscale*rand));
    v(1,t+1) = v(1,t+1) + 0.5 * dt * ((-v(1,t+1) - gsyn1 * s(2,t-dd+1) * (v(1,t+1) - vsyn) + iall + iapp1) * heav(tt-(tr(1)+trefrac)) + gel * ((v(2,t+1) - v(1,t+1)) + bspike * heav((tr(2)+dtr)-tt))+ normrnd(rnoise, signoise) * iappp(dt*t,iscale*rand) + iapp(dt*t, istim));
    v(2,t+1) = v(2,t+1) + 0.5 * dt * ((-v(2,t+1) - gsyn2 * s(1,t-dd+1) * (v(2,t+1) - vsyn) + iall + iapp2) * heav(tt-(tr(2)+trefrac)) + gel * ((v(1,t+1) - v(2,t+1)) + bspike * heav((tr(1) + dtr)-tt))); 
    % + normrnd(rnoise, signoise) * iappp(dt*t, iscale*rand));
    
    %tr(1) = 0;
    %tr(2) = 0;
    tt = tt + dt;
    s(1,t+1) = heav((tr(1)+dtr)-tt)*(1-s(1,t))/tauf - (1-s(1,t))/taus;
    s(2,t+1) = heav((tr(2)+dtr)-tt)*(1-s(2,t))/tauf - (1-s(2,t))/taus;
end

% plotting 
subplot(2,1,1)
h1 = plot(dt*t_grid,v(1,:));
%set(h,"linewidth", 2); 
%axis([0 max(t_grid) -95 60])
xlabel('Time [ms]'); ylabel('Voltage [mv]');
title('Excitory neuron (1)');

subplot(2,1,2)
h1 = plot(dt*t_grid,v(2,:));
%set(h,"linewidth", 2); 
%axis([0 max(t_grid) -95 60])
xlabel('Time [ms]'); ylabel('Voltage [mv]');
title('Inhibitory neuron (1)');