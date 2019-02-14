function predictplot(y, p)%, x, ylab, xlab)
% PREDICTPLOT Simple plot of prediction with intervals.
%
% predictplot(y)
% predictplot(y, q)

    switch nargin
        case 1
            p = [0.05, 0.50, 0.95];
        case 2
        if (length(p) ~= 3)
            error('predictplot: p argument must have length 3')
        end
        otherwise
            error('predictplot: incorrect number of input arguments');
    end

    x = 1:size(y)(1);
    quants = quantile(y, p, 2);  % TODO(brews): Check that this quantile method is compatible between Python, R, and MATLAB.

    c = [0.7, 0.7, 0.95];
    fill([x, fliplr(x)]', [quants(:, 3); flipud(quants(:, 1))], c, 'edgecolo', c), hold on
    plot(x, quants(:, 2))
end
