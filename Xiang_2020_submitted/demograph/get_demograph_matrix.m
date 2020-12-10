%% Get a demograph matrix for the cell signals in the cellList
% - Input
%   1. cellList: Oufti cellList
%   2. varargin:
%       - signal
%           Field name of the signal for construction of the demograph
%           e.g., 'signal0', 'signal1', 'signal2'
%           Default: 'signal0'
%       - comparator
%           Handle to a function that maps a cellStruct to an integer for
%           ranking.
%           e.g., lambda expression @(x) x.area
%           Default: @(x) x.length
%       - interp_factor
%           Height (number of rows) of the outptu demograph matrix. The
%           demograph is made of N columns, where N is the number of cells
%           in the cellList. This function extrapolates the signal values
%           to a fixed size that is represented by the interp_factor.
%           Default: 101
%       - pole_trim
%           Percentage of the signals to be excluded from the cell poles.
%           For example, if pole_trim = 0.05, 10% of the signal from the
%           cell poles would be excluded during the extrapolation.
%           Default: 0
%
% - Output
%   1. res: An (interp_factor) by N matrix, where the interp_factor is defaulted
%   to 101, and N is the number of cells in the cellList. The normalized
%   signal of each individual cell is represented as a column.
%   2. ranks: the ranking used for sorting the cells
%
% - Author
%   Yingjie Xiang, CJW Lab, Stanford University

function [res,ranks] = get_demograph_matrix(cellList,varargin)

% Default parameters
signal = 'signal0';
comparator = @(x) x.length;
interp_factor = 101;
pole_trim = 0;
normalize = true;

% Parser for optional input arguments
for ii = 1:2:length(varargin)
    switch lower(varargin{ii})
        case 'signal'
            signal = varargin{ii+1};
        case 'comparator'
            comparator = varargin{ii+1};
        case 'interp_factor'
            interp_factor = varargin{ii+1};
        case 'pole_trim'
            pole_trim = varargin{ii+1};
        case 'normalize'
            normalize= varargin{ii+1};
        otherwise
            error('Invalid optional arguments');
    end
end

% Contained for all signals
res = cell(length(cellList.meshData),1);
% Container for all ranking scores
ranks = cell(length(cellList.meshData),1);
% Go through each frame
for frame = 1:length(cellList.meshData)
    % Contained for the values of query for all cells in this frame
    vqs = zeros(interp_factor,length(cellList.meshData{frame}));
    % Corresponding ranks for all cells in this frame
    % Note the rank was flagged to -1
    rank = ones(length(cellList.meshData{frame}),1).*-1;
    % Go through each cell in this frame
    for cc = 1:length(cellList.meshData{frame})
        % Get extra data
        cellList.meshData{frame}{cc} = getextradata(cellList.meshData{frame}{cc});
        % Extract the cell signal
        v = cellList.meshData{frame}{cc}.(signal);
        % Normalized by the segment volume
        v = v ./ cellList.meshData{frame}{cc}.stepvolume;
        % Normalized the signal by the max value
        if normalize
            v = (v - min(v)) ./ (max(v) - min(v));
        end
        % Extraploation here
        x = linspace(0,1,length(v));
        xq = linspace(0+pole_trim,1-pole_trim,interp_factor);
        vqs(:,cc) = interp1(x,v,xq,'linear','extrap');
        % Validiate if a ranking score is possible
        try
            rank(cc) = comparator(cellList.meshData{frame}{cc});
        catch
            warning('Failed to rank Frame %d, Cell %d', frame,cc);
        end
    end
    res{frame} = vqs;
    ranks{frame} = rank;
end
% Collect results from all frames
res = horzcat(res{:});
% Collect ranks from all fames
ranks = vertcat(ranks{:});
% Check for cells with -1 rank, i.e. no rank score could be calculated
kill = ranks == -1;
% Remove the corresponding cells with -1 rank
res(:,kill) = [];
ranks(kill) = [];
% Sort the ranks
[ranks,six] = sort(ranks,'ascend');
% Sort the results correspondingly
res = res(:,six);
% Clear NaN
remove = any(isnan(res),1);
res(:,remove) = [];
ranks(remove) = [];
end