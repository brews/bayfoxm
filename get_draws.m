function [a, b, tau] = get_draws(seasonal_seatemp, foram)
% GET_DRAWS Fetch MCMC trace draws
%
% [a b t] = get_draws(seasonal_seatemp)
% [a b t] = get_draws(seasonal_seatemp, foram)

    if (nargin < 1)
        error('get_draws: not enough input arguments, 1 is required');
    end

    % First, ID which trace file we need.
    tracedir = './tracedumps/';  # Not sure this approach will always work.
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

    # Put together the column IDs we need and use it to pull correct col from 
    # `traces`.
    target_vars = {'a'; 'b'; 'tau'};
    n_vars = size(target_vars, 1);
    idxs = [];
    for i = 1:n_vars
        target_col = target_vars{i, 1};
        if strcmp(foram_half, 'hierarchical')
            target_col = [target_vars{i, 1} '__' foram];
        end
        target_idx = find(strcmp(header, target_col));
        if isempty(target_idx)
            error(['get_draws: cannot find parameters for foram: ' foram]);
        end
        idxs = [idxs; target_idx];
    end

    % Octave has no 'matsplit' so...
    a = traces(:, idxs(1));
    b = traces(:, idxs(2));
    tau = traces(:, idxs(3));
end
