
function [AF_SMT_data] = get_SMT_data()

AF_SMT_data = {};
addpath("./SMT_data/");

% (1) M9GlyCAA

dataset_num = 3;

RNAP_nameW = 'RNAP_WT_fracs_raw';
ribo_nameW = 'Ribosome_WT_fracs_raw';

RNAP_nameZ = 'RNAP_ftsZ_fracs_raw';
ribo_nameZ = 'Ribosome_ftsZ_fracs_raw';

RNAP_nameC = 'RNAP_oriC_fracs_raw';
ribo_nameC = 'Ribosome_oriC_fracs_raw';

[WT_data] = bin_data(RNAP_nameW, ribo_nameW, dataset_num);
[ftsZ_data] = bin_data(RNAP_nameZ, ribo_nameZ, dataset_num);
[oriC_data] = bin_data(RNAP_nameC, ribo_nameC, dataset_num);

AF_SMT_data{1}.WT = WT_data;
AF_SMT_data{1}.ftsZ = ftsZ_data;
AF_SMT_data{1}.oriC = oriC_data;

% (2) M9Gly

RNAP_nameW2 = 'M9gly_RNAP_WT_fracs_raw';
ribo_nameW2 = 'M9gly_Ribosome_WT_fracs_raw';

RNAP_nameC2 = 'M9gly_RNAP_oriC_fracs_raw';
ribo_nameC2 = 'M9gly_Ribosome_oriC_fracs_raw';

[WT_data2] = bin_data(RNAP_nameW2, ribo_nameW2, dataset_num);
[oriC_data2] = bin_data(RNAP_nameC2, ribo_nameC2, dataset_num);

AF_SMT_data{2}.WT = WT_data2;
AF_SMT_data{2}.oriC = oriC_data2;

% (3) M9Ala

RNAP_nameW3 = 'M9Ala_RNAP_WT_fracs_raw';
ribo_nameW3 = 'M9Ala_Ribosome_WT_fracs_raw';

RNAP_nameC3 = 'M9Ala_RNAP_oriC_fracs_raw';
ribo_nameC3 = 'M9Ala_Ribosome_oriC_fracs_raw';

[WT_data3] = bin_data(RNAP_nameW3, ribo_nameW3, dataset_num);
[oriC_data3] = bin_data(RNAP_nameC3, ribo_nameC3, dataset_num);

AF_SMT_data{3}.WT = WT_data3;
AF_SMT_data{3}.oriC = oriC_data3;

%temp = AF_SMT_data{1}.WT.RNAP.tableC;
%figure; plot(temp(:,1), temp(:,2), '.');

end


% ======================================================================

function [output] = bin_data(RNAP_name, ribo_name, n)

%(1) Load RNAP data

RNAP = {};
RNAP.data = load(RNAP_name);
RNAP.table = {};

for j = 1:n
    
    temp1 = RNAP.data.areas{j};
    temp2 = RNAP.data.fracs{j};
    RNAP.table{j} = [temp1' temp2'];
    
end

RNAP.tableC = [];

for j = 1:n

    RNAP.tableC = [RNAP.tableC; RNAP.table{j}];
    
end

% (2) Load ribo data

ribo = {};
ribo.data = load(ribo_name);
ribo.table = {};

for j = 1:n
    
    temp1 = ribo.data.areas{j};
    temp2 = ribo.data.fracs{j};
    ribo.table{j} = [temp1' temp2'];
    
end

ribo.tableC = [];

for j = 1:n

    ribo.tableC = [ribo.tableC; ribo.table{j}];
    
end


% ----------------------------------------------------------------------
% (0) Determine cell size bin
% ----------------------------------------------------------------------

bin_size = 0.5;         % cell area
bin_num = 12/bin_size;  % number of bins

data_range = NaN(bin_num, 2); % min and max for each bin

for b = 1:bin_num
    
    data_range(b,1) = bin_size *(b-1);
    data_range(b,2) = bin_size * b;
    
end

% ----------------------------------------------------------------------
% (i) RNAP: combine data from all table and partition by cell area bin
% ----------------------------------------------------------------------

for b = 1:bin_num
    
    bmin = data_range(b,1);
    bmax = data_range(b,2);

    temp = RNAP.tableC;
    temp( (temp(:,1) < bmin), : ) = [] ;
    temp( (temp(:,1) >= bmax), : ) = [] ;
        
    RNAP.bindata{b} = temp;
    
end

RNAP.stat = NaN(bin_num,3);  % num, mean, SD, SE

for b = 1:bin_num
    
    RNAP.stat(b,1) = size(RNAP.bindata{b}, 1);
    RNAP.stat(b,2) = nanmean(RNAP.bindata{b}(:,2));
    RNAP.stat(b,3) = nanstd(RNAP.bindata{b}(:,2));    
    RNAP.stat(b,4) = RNAP.stat(b,3)./sqrt(RNAP.stat(b,1));
    
end


% ----------------------------------------------------------------------
% (ii) Ribosome: combine data from all table and partition by cell area bin
% ----------------------------------------------------------------------

for b = 1:bin_num
    
    bmin = data_range(b,1);
    bmax = data_range(b,2);

    temp = ribo.tableC;
    temp( (temp(:,1) < bmin), : ) = [] ;
    temp( (temp(:,1) >= bmax), : ) = [] ;
        
    ribo.bindata{b} = temp;
    
end

ribo.stat = NaN(bin_num,3);  % num, mean, SD, SE

for b = 1:bin_num
    
    ribo.stat(b,1) = size(ribo.bindata{b}, 1);
    ribo.stat(b,2) = nanmean(ribo.bindata{b}(:,2));
    ribo.stat(b,3) = nanstd(ribo.bindata{b}(:,2));    
    ribo.stat(b,4) = ribo.stat(b,3)./sqrt(ribo.stat(b,1));
    
end

output = {};
output.RNAP = RNAP;
output.ribo = ribo;

end

% ======================================================================
