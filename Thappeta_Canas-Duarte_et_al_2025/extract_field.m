%% Extract data from a specific field from Oufti cellList
%
% -Input-
%   - cellList: Oufti cellList
%   - fieldname: a cell array of names to locate the field to be extracted
%       For example: {'A','B','C'}, this means the field to be extracted is
%       cellStruct.A.B.C. An arbitary 'depth' of the field is allowed.
%   - collapse (optional): if the result if a cell array, whether to
%       collapse the cell array into a single numberic array. If not
%       collapsed, each element in the cell array represents result from a
%       single cell.
%
% -Output-
%   - data: numeric array of all data in the specified field, or cell array
%   if the dimension of the data in the specified field does not all match
%
% -Author-
%   Yingjie Xiang, 2020-10-15, Jacobs-Wagner Lab, Stanford Universiy


function data = extract_field(cellList,fieldname,varargin)

collapse = false;
if nargin == 3
    collapse = varargin{1};
    assert(islogical(collapse),'3rd argument has to be logical: either true or false');
end

% Remove the frames with no cell detected
cellList.meshData(cellfun(@isempty,cellList.meshData)) = [];

data = cellfun(@(cells) cellfun(@(cellStruct) helper(cellStruct,fieldname),...
    cells, 'UniformOutput',0),cellList.meshData, 'UniformOutput',0);
data = [data{:}];

if all(cellfun(@numel,data) <= 1)
    data = vertcat(data{:});
elseif collapse && iscell(data)
    try
        data = vertcat(data{:});
    catch
        try
            data = horzcat(data{:})';
        catch
            warning('Cannot collapse the result, check array dimensions');
        end
    end
end

end

function res = helper(cellStruct,fieldname)
% The curr field is intially at the top level
curr = cellStruct;
% Keep search with a growing depth
ii = 1;
while ii <= length(fieldname)
    if isfield(curr,fieldname{ii})
        curr = curr.(fieldname{ii});
        ii = ii + 1;
    else
        break; % Did not find the field
    end
end
field_exist = ii == length(fieldname) + 1;
% Check if the field exists
% If not, return an empty array
if field_exist
    res = curr;
else
    res = [];
end
end