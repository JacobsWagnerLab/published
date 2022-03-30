
% Append cc_array for the Rsc statistics of mother cell (if exist)
% and daughtor cells

function [cc_array] = append_Rsc_mother_offspring(cc_array)

for m = 1:length(cc_array)

mother_ind = cc_array{m}.cc_info.mother_ind_on_cc_array;
offspringA_ind = cc_array{m}.cc_info.offspringA_ind_on_cc_array;
offspringB_ind = cc_array{m}.cc_info.offspringB_ind_on_cc_array;

Rsc_NaN = {};
Rsc_NaN.mean = NaN;
Rsc_NaN.std = NaN;
Rsc_NaN.max = NaN;
Rsc_NaN.min = NaN;
Rsc_NaN.maxdiff = NaN;

Rsc2_NaN = {};
Rsc2_NaN.trajR = NaN;
Rsc2_NaN.trajZ = NaN;
Rsc2_NaN.slope_R_cc_time = NaN;
Rsc2_NaN.slope_R_abs_time = NaN;
Rsc2_NaN.slope_Z_cc_time = NaN;
Rsc2_NaN.slope_Z_abs_time = NaN;

GR_NaN = {};
GR_NaN.GR = NaN;
GR_NaN.GR_CoD = NaN;
GR_NaN.ccT = NaN;
GR_NaN.midT = NaN;


if ( mother_ind > 0 )
    cc_array{m}.Rsc_statM = cc_array{mother_ind}.Rsc_stat;
    cc_array{m}.Rsc_statM2 = cc_array{mother_ind}.Rsc_stat2;
    cc_array{m}.GR_statM = cc_array{mother_ind}.GR_stat;
else
    cc_array{m}.Rsc_statM = Rsc_NaN;
    cc_array{m}.Rsc_statM2 = Rsc2_NaN;
    cc_array{m}.GR_statM = GR_NaN;
end

if ( offspringA_ind > 0)
    cc_array{m}.Rsc_statDA = cc_array{offspringA_ind}.Rsc_stat;
    cc_array{m}.Rsc_statDA2 = cc_array{offspringA_ind}.Rsc_stat2;
    cc_array{m}.GR_statDA = cc_array{offspringA_ind}.GR_stat;
else
    cc_array{m}.Rsc_statDA = Rsc_NaN;
    cc_array{m}.Rsc_statDA2 = Rsc2_NaN;
    cc_array{m}.GR_statDA = GR_NaN;
end

if ( offspringB_ind > 0)
    cc_array{m}.Rsc_statDB = cc_array{offspringB_ind}.Rsc_stat;
    cc_array{m}.Rsc_statDB2 = cc_array{offspringB_ind}.Rsc_stat2;
    cc_array{m}.GR_statDB = cc_array{offspringB_ind}.GR_stat;
else
    cc_array{m}.Rsc_statDB = Rsc_NaN;    
    cc_array{m}.Rsc_statDB2 = Rsc2_NaN;   
    cc_array{m}.GR_statDB = GR_NaN;  
end

end
