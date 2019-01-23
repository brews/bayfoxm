function ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram, drawsfun)
% PREDICT_SEATEMP Predict sea-surface temperature given δ18O of foram calcite and seawater δ18O.
%
% Syntax:
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std)
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp)
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram)
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram, drawsfun)
%
% Inputs:
%   d18oc - Scalar or [n x 1] column vector of δ18O of planktic foraminiferal 
%       calcite (‰; VPDB) samples.
%   d18osw - Scalar or n-length vector containing seatwater δ18O values (‰ VSMOW).
%   prior_mean - Scalar prior mean of sea-surface temperature (°C).
%   prior_std - Scalar prior standard deviation of sea-surface temperature (°C).
%   seasonal_seatemp - Boolean indicating whether to use annual or seasonal 
%       calibration model parameters.
%   foram - Character string. Foraminiferal of returned prediction. Can be 
%       "T. sacculifer", "N. pachyderma", "G. bulloides", "N. incompta", 
%       "G. ruber", or "none". If given species name, hierarchical calibration 
%       model is used. If not given or "none" uses pooled calibration model.
%   drawsfun - For testing, debugging and advanced use. Function used to get 
%       alpha, beta, tau Bayesian model parameters draw column vectors. Will be 
%       passed `seasonal_seatemp` and `foram` as arguments.
%
% Outputs:
%   ensemble - [n x m] matrix of sea-surface temperature (°C) predictions. 
%       m-length posterior ensemble for each of the n input sea-temperature 
%       values.
%
% Example:
%   % Annual pooled calibration model.
%   ens = predict_seatemp([0; 0.15; 0.2], 0.75, 15.0, 10.0)
%
%  % Hierarchical annual model for G. bulloides samples.
%   ens = predict_seatemp([0; 0.15; 0.2], 0.75, 15.0, 10.0, false, "G. bulloides")
%
%  % Hierarchical seasonal model for G. bulloides samples.
%   ens = predict_seatemp([0; 0.15; 0.2], 0.75, 15.0, 10.0, true, "G. bulloides")
%
% Other m-files required: get_draws.m
% Subfunctions: TARGET_ALL_PREDICT
% MAT-files required: None.
%
% See also: PREDICT_D18OC

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
    [a, b, tau] = drawsfun(seasonal_seatemp, foram);
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
% TARGET_ALL_PREDICT Analytical solution for seatemp inference with single MCMC draw.
%
% Syntax: y = target_all_predict(d18oc, d18osw, a, b, tau2, prior_mu, prior_inv_cov)
%
% Inputs:
%   d18oc - Scalar or [n x 1] column vector of δ18O of planktic foraminiferal 
%       calcite (‰; VPDB) samples.
%   d18osw - Scalar or n-length column vector containing seatwater δ18O values (‰ VPDB).
%   a - Scalar MCMC trace draw of intercept parameter.
%   b - Scalar MCMC trace draw of slope parameter.
%   tau2 - Scalar MCMC trace draw of residual variance.
%   prior_mu - [n, 1] column vector of prior means for each element of `d18oc`.
%   prior_inv_cov - [n, n] matrix, inverse of the prior covariance matrix for 
%       `d18oc`.
%
% Outputs:
%   y - Vector of sample inferred seatemp time series conditional on input.

    % NOTE: δ18Osw should already be scale adjusted to VPDB from VSMOW.
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
