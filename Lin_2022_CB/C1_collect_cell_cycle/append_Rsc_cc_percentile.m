
% Normalize cell cycle and examine the Rsc dynamics

function [cc_array] = append_Rsc_cc_percentile(cc_array)

for m = 1:length(cc_array)

cc_array{m}.Rsc_stat2.trajR = NaN;
cc_array{m}.Rsc_stat2.trajZ = NaN;
cc_array{m}.Rsc_stat2.slope_R_cc_time = NaN;
cc_array{m}.Rsc_stat2.slope_R_abs_time = NaN;
cc_array{m}.Rsc_stat2.slope_Z_cc_time = NaN;
cc_array{m}.Rsc_stat2.slope_Z_abs_time = NaN;

% Rsc trajectory and the measurement time points
time = cc_array{m}.signal(:,1);
Rsc = cc_array{m}.Rsc_traj;

if ( length(Rsc) > 2 )
    
% initial and ending time of cell cycle
ini_time = cc_array{m}.data(1,1);
end_time = cc_array{m}.data(end,1);

% convert into cell cycle percentile
time_cc_percent = (time - ini_time) / (end_time-ini_time);
t = ((1:100)/100)';
Rsc_cc_fit = interp1(time_cc_percent, Rsc, t, 'pchip');

% linear fitting Rsc to get slope 
coefR = polyfit(t, Rsc_cc_fit, 1);
slope_R_cc_time = coefR(1);
slope_R_abs_time = coefR(1)/(end_time-ini_time);

% z-transformed Rsc
Rsc_cc_ZT = (Rsc_cc_fit - mean(Rsc_cc_fit))/ std(Rsc_cc_fit);
coefZ = polyfit(t, Rsc_cc_ZT, 1);
slope_Z_cc_time = coefZ(1);
slope_Z_abs_time = coefZ(1)/(end_time-ini_time); 

% Save the result

cc_array{m}.Rsc_stat2.trajR = Rsc_cc_fit;
cc_array{m}.Rsc_stat2.trajZ = Rsc_cc_ZT;
cc_array{m}.Rsc_stat2.slope_R_cc_time = slope_R_cc_time;
cc_array{m}.Rsc_stat2.slope_R_abs_time = slope_R_abs_time;
cc_array{m}.Rsc_stat2.slope_Z_cc_time = slope_Z_cc_time;
cc_array{m}.Rsc_stat2.slope_Z_abs_time = slope_Z_abs_time;

end

end

%figure;
%plot(time_cc_percent, Rsc, 'o'); hold on;
%plot(t, Rsc_cc_fit, '.-');

end
