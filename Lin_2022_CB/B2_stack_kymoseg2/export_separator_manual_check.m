
function [sep_summary] = export_separator_manual_check(separatorC, export_dir)

% Export an excel sheet for summary of separatorC 

L = length(separatorC);

sep_summary = zeros(L,5);
% (1) first frame, (2) last frame (3) total span of frame
% (4) separator check by default  (5) manual check

for s = 1:L

    temp = separatorC{s};

    ini_frame = temp.data(1,1);
    end_frame = temp.data(end,1);
    total_frame = size(temp.data,1);    
    sep_check = temp.length && (temp.ending || temp.exit_y);
    
    write = [ini_frame end_frame total_frame sep_check sep_check];

    sep_summary(s,:) = write;    
    
end

cd(export_dir);
writematrix(sep_summary, 'sep_manual_check.xls', 'WriteMode', 'overwritesheet');

