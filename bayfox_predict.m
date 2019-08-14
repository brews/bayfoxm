function output = bayfox_predict(d18oc,d18osw,prior_mean,prior_std,species,bayes)
% function output = bayfox_predict(d18oc,d18osw,prior_mean,prior_std,species,bayes)
%
% BAYFOX prediction model for d18O of planktic foraminifera.
% Predicts SST from d18O of foraminiferal calcite, with given d18O of
% seawater, and species group.
% ----- Inputs -----
% d18oc: A scalar or vector of d18O of foram calcite (in VPDB) (1 x N) or (N x 1)
%
% d18osw: A scalar or vector of d18O of seawater (in VSMOW) (1 x N) or (N x 1)
%
% prior_mean: A scalar prior mean value of SST in degrees C.
%
% prior_std: A scalar prior standard deviation value of SST in degrees C.
%
% species: A character string corresponding the foram species. Options are:
% 'bulloides' = G. bulloides
% 'incompta' = N. incompta
% 'pachy' = N. pachyderma
% 'ruber' = G. ruber
% 'sacculifer' = T. sacculifer
% 'all' = use the pooled annual (non-species specific) model
% 'all_sea' = use the pooled seasonal (non-species specific) model
%
% bayes: (optional - mainly for use with the DASH interface) A Bayesian
% posterior to use for the calibration. if empty, the default posterior
% file is used.
%
% ----- Outputs -----
%
% output is a structure containing:
% output.prior_mean: User choice of prior mean
% output.prior_std: User choice of prior standard deviation
% output.SST: 2.5%, 50%, and 97.5% confidence intervals of posterior SST
% output.ens: full ensemble of posterior SST (N x 2000);
%
    % Ensure column vectors.
    d18oc=d18oc(:);
    d18osw=d18osw(:);
    species_list = {'bulloides','incompta','pachy','ruber','sacculifer','all','all_sea'};
    %check that you have an allowable species
    if ~ismember(species,species_list)
        error('Species not recognized');
    end
    % process optional call to a different posterior file
    ng=nargin;
    if ng == 6
    else
        bayes=["poolann_params.mat";"poolsea_params.mat";"hiersea_params.mat"];
    end
    %
    id = (1:1:5);
    %load appropriate model
    if  strcmp(species,'all')
        params = load(bayes(1));
        id = 1;
    elseif strcmp(species,'all_sea')
        params = load(bayes(2));
        id = 1;
    else
        params = load(bayes(3));
        %grab id location for species.
        id = id(ismember(species_list,species));
    end
    %get posterior params
    betaT=params.betaT(:,id);
    alpha=params.alpha(:,id);
    sigma=params.sigma(:,id);
    %get dimensions of time series and draws
    nd = length(d18oc);
    n_draws = length(sigma);
    % Unit adjustment for permil VSMOW to permil VPDB.
    d18osw_adj = d18osw - 0.27;
    % Prior mean and inverse covariance matrix
    pmu = repmat(ones(nd, 1) * prior_mean,1,n_draws);
    pinv_cov = repmat(prior_std,nd,n_draws).^-2;
    
    % Posterior calculations
    post_mean_num = pinv_cov .* pmu + repmat(sigma',nd,1).^-2 .* repmat(betaT',nd,1) .* (d18oc - repmat(alpha',nd,1) - d18osw_adj);
    post_mean_den = pinv_cov + repmat(betaT',nd,1).^2 .* repmat(sigma',nd,1).^-2;
    post_mean = post_mean_num ./ post_mean_den;
    post_sig = sqrt(post_mean_den.^-1);
    output.ens = post_mean + randn(nd,n_draws).*post_sig;
    output.prior_mean = prior_mean;
    output.prior_std = prior_std;
    % truncate at T < -2.5 (seawater is frozen)
    output.ens(output.ens < -2.5) = NaN;
    output.SST = prctile(sort(output.ens,2),[2.5 50 97.5],2);
end