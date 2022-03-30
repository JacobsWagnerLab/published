
function [mother_cell] = append_Rsc_sensor(mother_cell, param)

Rsc_adjust_factor = param.Rsc_adjust_factor;
sensor_adjust_factor = param.sensor_adjust_factor;
exposure_adjust_factor = param.exposure_adjust_factor;

% (1) Coonstruct Rsc and sensor trajectory 
%    (optional) Remove the data after last division

Rsc = {};
sensor = {};

%mother_cell.vec( isnan(mother_cell.vec(:,1)), : ) = [];

last_div_cc = mother_cell.vec(end,1);
last_div_frame = find(mother_cell.vec(:,1)== last_div_cc, 1, 'first');

% [1] time [2] c2 channel [3] c3 channel
time = (1:length(mother_cell.signal))';
signal = [time mother_cell.signal];
signal( isnan(signal(:,2)), : ) = [];

Rsc.data = [signal(:,1) Rsc_adjust_factor* ( signal(:,2)./signal(:,3) )];
sensor.data = [signal(:,1) sensor_adjust_factor(1)*signal(:,2) + sensor_adjust_factor(2)*signal(:,3)];
sensor.data(:,2) = sensor.data(:,2)/ exposure_adjust_factor ;

% Truncate the trajectory at last cell division
Rsc.dataS = Rsc.data;
Rsc.dataS( ( Rsc.dataS(:,1) > last_div_frame-1 ),: ) = [];
        
sensor.dataS = sensor.data;
sensor.dataS( ( sensor.dataS(:,1) > last_div_frame-1 ),: ) = [];

mother_cell.Rsc = Rsc;
mother_cell.sensor = sensor;

end


