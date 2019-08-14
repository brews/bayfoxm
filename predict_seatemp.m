function ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram)
% PREDICT_SEATEMP Predict sea-surface temperature given d18O of foram calcite and seawater d18O.
%
% Syntax:
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std)
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp)
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram)
%   ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std, seasonal_seatemp, foram, drawsfun)
%
% Inputs:
%   d18oc - Scalar or [n x 1] column vector of d18O of planktic foraminiferal 
%       calcite (permil VPDB) samples.
%   d18osw - Scalar or n-length vector containing seatwater d18O values (permil VSMOW).
%   prior_mean - Scalar prior mean of sea-surface temperature (°C).
%   prior_std - Scalar prior standard deviation of sea-surface temperature (°C).
%   seasonal_seatemp - Boolean indicating whether to use annual or seasonal 
%       calibration model parameters.
%   foram - Character string. Foraminiferal of returned prediction. Can be 
%       "T. sacculifer", "N. pachyderma", "G. bulloides", "N. incompta", 
%       "G. ruber", or "none". If given species name, hierarchical calibration 
%       model is used. If not given or "none" uses pooled calibration model.
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
% Subfunctions: None.
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
    %Note that tau is standard deviation
    
    %get dimensions of time series and draws
    nd = length(d18oc);
    n_draws = length(tau);

    % Unit adjustment for permil VSMOW to permil VPDB.
    d18osw_adj = d18osw - 0.27;

    % Prior mean and inverse covariance matrix
    pmu = repmat(ones(nd, 1) * prior_mean,1,n_draws);
    pinv_cov = repmat(prior_std,nd,n_draws).^-2;
    
    % Posterior calculations
    post_mean_num = pinv_cov .* pmu + repmat(tau',nd,1).^-2 .* repmat(b',nd,1) .* (d18oc - repmat(a',nd,1) - d18osw_adj);
    post_mean_den = pinv_cov + repmat(b',nd,1).^2 .* repmat(tau',nd,1).^-2;
    post_mean = post_mean_num ./ post_mean_den;
    post_sig = sqrt(post_mean_den.^-1);
    ensemble = post_mean + randn(nd,n_draws).*post_sig;
end