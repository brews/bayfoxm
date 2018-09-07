function ensemble = predict_d18oc(seatemp, d18osw, seasonal_seatemp, foram, drawsfun)
% PREDICT_D18OC Predict sea-surface temperature given δ18O of foram calcite and seawater δ18O.
%
% [alpha beta tau] = predict_d18oc(seatemp, d18osw)
% [alpha beta tau] = predict_d18oc(seatemp, d18osw, seasonal_seatemp)
% [alpha beta tau] = predict_d18oc(seatemp, d18osw, seasonal_seatemp, foram)
% [alpha beta tau] = predict_d18oc(seatemp, d18osw, seasonal_seatemp, foram, drawsfun)

    switch nargin  % Set default args. A bit janky.
        case 2
            seasonal_seatemp = false;
            foram = 'none';
            drawsfun = @get_draws;
        case 3
            foram = 'none';
            drawsfun = @get_draws;
        case 4
            drawsfun = @get_draws;
        case 5
            # DEBUG RUN
        otherwise
            error('predict_d18oc: incorrect number of input arguments');
    end

    % Get MCMC trace draws
    % This decides whether we use pooled vs seasonal and pooled vs hierarchical.
    [a, b, tau] = get_draws(seasonal_seatemp, foram);
    tau2 = tau .^ 2;

    nd = length(seatemp);
    n_draws = length(tau);

    % Unit adjustment for permille VSMOW to permille VPDB.
    d18osw_adj = d18osw - 0.27;

    # Might be able to vectorize the loop below...
    ensemble = NaN(nd, n_draws);
    for (i = 1:n_draws)
        ensemble(:, i) = normrnd(a(i) + seatemp * b(i) + d18osw_adj, tau2(i));
    end
end
