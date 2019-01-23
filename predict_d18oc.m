function ensemble = predict_d18oc(seatemp, d18osw, seasonal_seatemp, foram, drawsfun)
% PREDICT_D18OC - Predict δ18O of foram calcite given sea-surface temperature and seawater δ18O.
%
% Syntax:
%   ensemble = predict_d18oc(seatemp, d18osw)
%   ensemble = predict_d18oc(seatemp, d18osw, seasonal_seatemp)
%   ensemble = predict_d18oc(seatemp, d18osw, seasonal_seatemp, foram)
%   ensemble = predict_d18oc(seatemp, d18osw, seasonal_seatemp, foram, drawsfun)
%
% Inputs:
%   seatemp - Scalar or [n x 1] column vector of sea-surface temperature values (°C).
%   d18osw - Scalar or n-length vector containing seatwater d18O values (‰ VSMOW).
%   seasonal_seatemp - Boolean indicating whether to use annual or seasonal 
%       calibration model parameters.
%   foram - Character string. Foraminiferal of returned prediction. Can be 
%       "T. sacculifer", "N. pachyderma", "G. bulloides", "N. incompta", 
%       "G. ruber", or "none". If given species name, hierarchical calibration 
%       model is used.  If not given or "none" uses pooled calibration model.
%   drawsfun - For testing, debugging and advanced use. Function used to get 
%       alpha, beta, tau Bayesian model parameters draw column vectors. Will be 
%       passed `seasonal_seatemp` and `foram` as arguments.
%
% Outputs:
%   ensemble - [n x m] matrix of δ18Oc (‰ VPDB) predictions. m-length posterior 
%       ensemble for each of the n input sea-temperature values.
%
% Example:
%   % Annual pooled calibration model.
%   ens = predict_d18oc([15; 15.5; 16], 0.75)
%
%   % Hierarchical annual model for G. bulloides samples.
%   ens = predict_d18oc([15; 15.6; 16], 0.75, false, "G. bulloides")
%
%   % Hierarchical seasonal model for G. bulloides samples.
%   ens = predict_d18oc([15; 15.6; 16], 0.75, true, "G. bulloides")
%
% Other m-files required: get_draws.m
% Subfunctions: None.
% MAT-files required: None.
%
% See also: PREDICT_SEATEMP

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
