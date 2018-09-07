function ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram, drawsfun)
% PREDICT_D18OC Predict sea-surface temperature given δ18O of foram calcite and seawater δ18O.
%
% [alpha beta tau] = predict_d18oc(d18oc, d18osw, prior_mean, prior_std)
% [alpha beta tau] = predict_d18oc(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp)
% [alpha beta tau] = predict_d18oc(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram)
% [alpha beta tau] = predict_d18oc(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram, drawsfun)

    switch nargin  % Set default args. A bit janky.
        case 4
            seasonal_seatemp = false;
            foram = 'none';
            drawsfun = @get_draws;
        case 5
            foram = 'none';
            drawsfun = @get_draws;
        case 6
            drawsfun = @get_draws;
        otherwise
            error('predict_seatemp: incorrect number of input arguments');
    end

    % Get MCMC trace draws
    % This decides whether we use pooled vs seasonal and pooled vs hierarchical.
    [a, b, tau] = get_draws(seasonal_seatemp, foram);
    tau2 = tau .^ 2;

    nd = length(d18oc);
    n_draws = length(tau);

    % Unit adjustment for permille VSMOW to permille VPDB.
    d18osw_adj = d18osw - 0.27;

    % Prior mean and inverse covariance matrix
    pmu = ones(nd, 1) * prior_mean;
    pinv_cov = eye(nd) * prior_std ^ (-2);

    # Might be able to vectorize the loop below...
    ensemble = NaN(nd, n_draws);
    for (i = 1:n_draws)
        ensemble(:, i) = target_all_predict(d18oc, d18osw, a(i), b(i), tau2(i), pmu, pinv_cov);
    end
end


function y = target_all_predict(d18oc, d18osw, a, b, tau2, prior_mu, prior_inv_cov)
% TARGET_ALL_PREDICT Analytical solution for d18Oc inference

    % d18osw should be scale adjusted to VPDB from VSMOW.
    if (nargin ~= 7)
        error('target_all_predict: incorrect number of input arguments');
    end

    n_ts = length(d18oc);

    % Inverse posterior covariance matrix
    inv_post_cov = prior_inv_cov + b ^ 2 / tau2 * eye(n_ts);

    % Inverse to get posterior covariance matrix, following BAYSPAR, we'll use
    % cholesky
    opts.SYM = true;
    opts.POSDEF = true;
    post_cov = linsolve(inv_post_cov, eye(n_ts), opts);

    % Get square root via cholesky factor
    sqrt_post_cov = chol(post_cov)';

    % mean first factor
    mean_ff = prior_inv_cov * prior_mu + (1 / tau2) * b * (d18oc - a - d18osw);
    mean_full = post_cov * mean_ff;

    y = mean_full + sqrt_post_cov * randn(n_ts, 1);
end
