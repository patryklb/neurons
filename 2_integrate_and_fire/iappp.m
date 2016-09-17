function res = iappp(t, pon)
    load('params.mat')
    res = heav(floor(t)+d-t).* heav(t-floor(t)) .* heav(t-pon);
end
