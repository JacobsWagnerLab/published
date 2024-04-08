
% =======================================================================
% Extract cell cycle parameters from cell size distributions of each type
% % =====================================================================

path = {};
folder = '(dir)./save_20230509B/';
path{1} = strcat(folder, '0405/');
path{2} = strcat(folder, '0410/');
path{3} = strcat(folder, '0412/');
file_name = 'Type_summary.mat';

% (1) Analyze each dataset separately

GR_data = {};
data = {};

for DB = 1:3
    
    cd(path{DB});
    data{DB} = load(file_name).rec_data;

    rec = sub_analyze_CCG_V2(data{DB});
    GR_data{DB} =  GR_analysis_CCG_V2(rec);

end

%

rec = {};

for DB = 1:3
    
    rec.Area(DB,:) = GR_data{DB}.Area;
    rec.Count(DB,:) = GR_data{DB}.type_count;
    rec.Mnull(DB,:) = GR_data{DB}.Mnull;
    rec.M(DB,:) = GR_data{DB}.M;
    rec.muL(DB,:) = GR_data{DB}.muL;
    
end
    

% (2) Pool all data and analyze as single dataset

data_combine = [];

for DB = 1:3    
    data_combine = [data_combine; data{DB}];
end

rec_combine = sub_analyze_CCG_V2(data_combine);
GR_data_combine = sub_GR_analysis_CCG_V2(rec_combine);

