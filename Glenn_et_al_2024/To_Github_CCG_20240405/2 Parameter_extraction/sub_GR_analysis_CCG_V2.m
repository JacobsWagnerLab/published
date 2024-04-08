% ========================================================================
% Extract realtive fraction, relative growth rate, relative cell cycle time
% of each cell cycle stages
% ========================================================================

function GR_data = sub_GR_analysis_CCG_V2(rec)

type_count = rec.type_count;
Area = rec.Area;

total_num = sum(type_count);
M(1) = type_count(1) /total_num;
M(2) = type_count(2) /total_num;
M(3) = type_count(3) /total_num;

Mnull(1) = 1 - (Area(1)/Area(2));
Mnull(2) = ( (Area(3)/Area(2)) - 1 ) *  ( (Area(1)+Area(2))/Area(3) );
Mnull(3) = 1 - Mnull(1) - Mnull(2);

[muL, TL] = Sub_Three_state_model_CCG_V2(Area, M);

GR_data = {};
GR_data.Area = Area;
GR_data.type_count = type_count;

GR_data.M = M;
GR_data.Mnull = Mnull;
GR_data.muL = muL';
GR_data.TL = TL';

% ==================================================================

% Use 3-state model to calculate stage-specific elongation rate

function [muL, TL] = Sub_Three_state_model_CCG_V2(Area, M)

%M = [0.4846 0.2671 0.2482];
%Area = [250 300 400];

TL = NaN(3,1);  % tau* lambda
TL(1) = log( 1 / (1-M(1)) );
TL(2) = log( (1 +M(2) + M(3)) / ( 1+M(3) ) );
TL(3) = log( 1 + M(3) );

LR = NaN(3,1);  % log ratio
LR(1) = log( Area(2)/Area(1) );
LR(2) = log( Area(3)/Area(2) );
LR(3) = log( (Area(1)+Area(2)) /Area(3) );

muL = LR./TL;

% ==================================================================
