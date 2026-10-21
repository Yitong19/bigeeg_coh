%% resampling function (monotone Hermite spline)
function c_new = interpolate_column(t, c, t_new)

    % Keep only finite time/value pairs
    ok = isfinite(t) & isfinite(c);
    t_valid = t(ok);
    c_valid = c(ok);

    % Sort timestamps and remove duplicates, keeping the first value
    [t_valid, idx] = unique(t_valid, 'first');
    c_valid = c_valid(idx);

    % Shape-preserving cubic interpolation
    c_new = interp1(t_valid, c_valid, t_new, 'pchip', 'extrap');
end