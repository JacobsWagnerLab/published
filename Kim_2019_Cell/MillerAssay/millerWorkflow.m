% One should have Excel output from a microplate reader
% First column = time points of reading (in minutes)
% second column = OD420 reading in the first well in each time point.
rawData = 0;
% copy and paste the Excel output into rawData in Workspace.
% Usually the first column is...
% timeseries = [9;19;29;39;49;59;69;79;89;99;109;119;129;149;169;189;209;229;249;269;289;309;329;349;369;389;409;429;449;469;489;509;529;549;569;589;609;629;649;669;689;709;729;749;769;789;809;829;849;869;];
mslope = MillerslopeKH(rawData,2);


% Alternative approach in calculating mslope;
% rawData, same as above
od550 = 0;
% For od550, copy and paste Excel output into od550, 
% it is OD550 reading in each well for each timepoint when OD420 was
% measured.
mslope = MillerslopeWolf(rawData,od550,2);

% for LacZ activity, I used this equation. Calculation was done in Excel.
% MU = 1000*mslope/OD600/0.597;

% Once I have MU vs time, sqrt(MU) vs time, fit with y = c and y = a+bx as
% written in the Kim et al paper (materials and methods)
% from each day, I will have (a, b, c). save them in X.
% for 0.05 mM raw 
% X = [1.24680000000000,0.228000000000000,1.63460000000000;1.07520000000000,0.311200000000000,1.60880000000000;-0.685300000000000,0.895890000000000,1.17390000000000;0.823410000000000,0.268190000000000,1.30660000000000;-0.664700000000000,1.10530000000000,1.17100000000000;0.381070000000000,0.542980000000000,1.40620000000000;0.324090000000000,0.607020000000000,1.55190000000000;1.11130000000000,0.206720000000000,1.48540000000000;0.658200000000000,0.382390000000000,1.52350000000000;0.752770000000000,0.241390000000000,1.10870000000000;0.595420000000000,0.469420000000000,1.33360000000000;0.662650000000000,0.421980000000000,1.34780000000000;0.924080000000000,0.174270000000000,1.26630000000000];
fit = Millerfitavg(X);

