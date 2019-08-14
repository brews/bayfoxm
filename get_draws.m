function [a, b, tau] = get_draws(seasonal_seatemp, foram)
% GET_DRAWS - Fetch MCMC trace draws
% This is used internally by `predict_d18oc` and `predict_seatemp`.
%
% Syntax:
%   [a b tau] = get_draws(seasonal_seatemp)
%   [a b tau] = get_draws(seasonal_seatemp, foram)
%
% Inputs:
%   seasonal_seatemp - Boolean. Indicates whether to return draws from a 
%       seasonal model. If `false`, returns annual model draws. Default is 
%       `false.
%   foram - Optional. Character string. Which foram group to return draws for, 
%       in a hierarchical model. If not given or "none", will use pooled model.
%
%  Outputs:
%   a - [n x 1] column vector of alpha or intercept parameter for Bayesian model.
%   b - [n x 1] column vector of beta or slope parameter for Bayesian model.
%   tau - [n x 1] column vector of tau or error term of Bayesian model.
%
% Example:
%   [a b tau] = get_draws(false)  % Annual pooled model parameters.
%
%   % Hierarchical annual model for G. bulloides samples.
%   [a b tau] = get_draws(false, "G. bulloides")
%
%   % Hierarchical seasonal model for G. bulloides samples.
%   [a b tau] = get_draws(true, "G. bulloides")
%
% Other m-files required: None.
% Subfunctions: None.
% MAT-files required: None.
%
% Depends on ./tracedumps/*.csv directory and csv file names. Be sure that this 
% directory is in PWD.

    if (nargin < 1)
        error('get_draws: not enough input arguments, 1 is required');
    end

    % Get path to where the trace files are
    path_ind=which('get_draws');
    path_ind=path_ind(1:end-11);
    tracedir = strcat(path_ind,'tracedumps/');
    % ID which trace file we need.
    foram_half = 'hierarchical';
    time_half = 'annual';

    if (nargin == 1 || isempty(foram) || strcmp('none', foram))
        foram_half = 'pooled';
    end

    if (seasonal_seatemp == true)
        time_half = 'seasonal';
    end

    tracepath = [tracedir time_half '_' foram_half '_trace.csv'];

    % Read file traces
    traces = csvread(tracepath, 1, 0);
    % Now get file headers
    fid = fopen(tracepath);
    first_line = fgetl(fid);
    fclose(fid);
    header = strsplit(first_line, ',');

    %Put together the column IDs we need and use it to pull correct col from 
    % `traces`.
    target_vars = {'a'; 'b'; 'tau'};
    n_vars = size(target_vars, 1);
    idxs = NaN(3,1);
    for i = 1:n_vars
        target_col = target_vars{i, 1};
        if strcmp(foram_half, 'hierarchical')
            target_col = strcat(target_vars{i, 1},'__',foram);
        end
        target_idx = find(header==target_col);
        if isempty(target_idx)
            error(['get_draws: cannot find parameters for foram: ' foram]);
        end
        idxs(i) = target_idx;
    end

    % Octave has no 'matsplit' so...
    a = traces(:, idxs(1));
    b = traces(:, idxs(2));
    tau = traces(:, idxs(3));
    % thin these draws to 1000
    a = a(1:10:end);
    b = b(1:10:end);
    tau = tau(1:10:end);
end
