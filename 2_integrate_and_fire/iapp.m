function res = iapp(t, istim)
    load('params.mat')
    res = heav(pon+deltat-t).*heav(t-pon) * istim;
end

