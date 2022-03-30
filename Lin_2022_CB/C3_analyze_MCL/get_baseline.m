
function [baseline] = get_baseline(time, data, method_type, param)

% Getting baseline of a trajectory uing different methods

baseline = [];

if (method_type == 1)

    % Method 1: Using local convolution with specified window
    window_size = param;
    window = ones(window_size, 1) / window_size;

    % Using reflexive boundary condition to get convolution at boundary
    data_expand = [flip(data); data; flip(data)];
    baseline = conv(data_expand, window, 'same');

    L = length(time);
    baseline = baseline(L+1:2*L);

elseif (method_type == 2)

    % Method 2: Using polynomial fit

    poly_degree = param;

    coef = polyfit(time, data, poly_degree);
    baseline = polyval(coef, time);

end

