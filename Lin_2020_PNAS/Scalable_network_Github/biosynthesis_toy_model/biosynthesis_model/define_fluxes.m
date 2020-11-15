
% This script describes the formula of flux functions (J), and the
% parameter information (param). These information are returned to the 
% ODE model for integrating the system. 

% Each flux function has a formula (in string format) and associated parameter(s). 

% The formula is described in 'J'. 
% There are 40 flux fucntion in this model, hence 'J' is a cell array containing 24 strings. 

% The information a given parameter is described in 'param'. 
% There are total 78 parameters in this model, hence 'param' is a 78-by-3 matrix, 
% each row specifies information of one flux. 

% The columns specifies following information:
% column [1]: Index of flux function being associated
% column [2]: The rank (first, second, etc) of this parameter in 
%             the associated flux function.
% column [3]: The type of the associated flux function.


function [J, output] = define_fluxes()

% There are 40 flux functions in this model, belong to five "TYPES". 
% Flux functions with the same "TYPES" have same algebraic structure, 
% the difference is their variable are different. 

% For example, for the Constant Degradation (CD) TYPE, 
% it is described by a linear function q(i)*u(j)
% where q(i) is the i'th parameter, 
% and u(j) is the j'th variable (corresponded to j'th node in the system)

% There are five TYPES of algebraic structure in this model: 
% (1) CD, (2) MS1, (3) MS2, (4) MS3, (5) Mpoly (see description in below)
% the formula of these TYPES are based on Michaelis-Menten kinetics 
% and is described in Supplementary Materials. 

% Here, instead of entering formula for 24 flux function, 
% we use "TYPE" as a template to generate formula for all flux function. 

% Each "TYPE" below is a template string, eg. 'q(1)*u(1)'
% When defining a function, this string is replaced by corresponded
% variable and parameters to generate a formula string of a flux function.

% For example, for Constant Degredation flux function with parameter p(i) and variable y(j),
% the string 'q(1)*u(1)' is replaced by 'p(i)*y(j)'. This step is performed
% by subfunction 'ind_assignV.m'. 


% ====== Define TYPES for flux function ==================================

% Belows are all string operation to construct the formula of flux functions

% (i) MS1: single-substrate MMS
MS1 = 'q(1)*u(1)* (u(2)/( u(2)+ q(2)) )';

% (ii) MS2: single-substrate MMS with ATP
MS2 = 'q(1)*u(1)*u(2)*u(3) / (q(3)*q(2) + q(2)*u(2) + u(2)*u(3) ) ';

% (iii) MS3: bi-substrate MMS with ATP
MS3 = 'q(1)*u(1)*u(2)*u(3)*u(4) / (q(4)*q(3)*q(2) + q(4)*q(2)*u(2) + q(2)*u(2)*u(3) + u(2)*u(3)*u(4) )';

% (iv) Poly: polymerization course-grained reaction

MP1 = '( q(3)*q(4) + q(3)*u(4) + u(2)*u(4) )';
MP2 = '(u(4)/u(5)) *( q(3)*q(4) + q(3)*u(5) + u(2)*u(5) )';
MP3 = '(u(4)/u(6)) *( q(3)*q(4) + q(3)*u(6) + u(2)*u(6) )';
MP4 = '(u(4)/u(7)) *( q(3)*q(4) + q(3)*u(7) + u(2)*u(7) )';

Mdeo = strcat('(q(1)/u(3))*', MP1, '+ q(5)*', MP1, '+q(6)*', MP2, '+q(7)*', MP3, '+q(8)*', MP4 );
Mnum = 'q(2)*u(1)*u(2)*u(4)';
Mpoly = strcat(Mnum, '/(', Mdeo, ')');

% (v) CD: constant degradation 
CD = 'q(1)*u(1)';

% ============================================================== %

% Quick reference for node index: 

% [1-10]  m1-m10 (metabolites)
% [11]    ATP
% [12]    ADP
% [13-15] P1-P3  (transporters)
% [16-22] Q1-Q7  (house-keeping enzymes)
% [23]    R      (ribosome)
% [24]    Q8     (template such as RNA)

%%% ================ Specifies flux function ================= %%%

global pcount;
pcount = 0;   % current parameter index

global Jcount;
Jcount = 1;   % corrent flux function index

global param; 
param = [];

J = {};

% Quick reference of flux TYPE 
% (i) MS1, (ii) MS2, (iii) MS3, (iv) Poly (v) CD

% (1) Ain-Cin (type MS1 = 1)
% These are transporter fluxes that depends on transporters (#13-#15)
% and ATP as energy input (#11)

J{Jcount} = ind_assignV( MS1, [13 11], 1 );
J{Jcount} = ind_assignV( MS1, [14 11], 1 );
J{Jcount} = ind_assignV( MS1, [15 11], 1 );

% (2) S1-S3 (type MS3 = 3)
% These are binary synthesis fluxes that depends on enzymes (#16-18),
% metabolites (#1-3) and ATP (#11)

J{Jcount} = ind_assignV( MS3, [16 1 2 11], 3 );
J{Jcount} = ind_assignV( MS3, [17 1 3 11], 3 );
J{Jcount} = ind_assignV( MS3, [18 2 3 11], 3 );

% (3a) R1F, R2F (type MS3 = 3)
% These are reversible forward-conversion fluxes, that depends on the
% enzymes (#19-20), metabolites (#4-6) and ATP (#11)

J{Jcount} = ind_assignV( MS3, [19 4 5 11], 3 );
J{Jcount} = ind_assignV( MS3, [20 5 6 11], 3 );

% (3b) Define R1R, R2R (type MS2 = 2)
% These are reversible reverse -conversion fluxes, that depends on the
% enzymes (#19-20), metabolites (#7-10)

J{Jcount} = ind_assignV( MS2, [19 7 8], 2 );
J{Jcount} = ind_assignV( MS2, [20 9 10], 2 );

% (4a) SATP (type MS1 = 1) 
% This is the de novo ADP-synthesis flux, depends on enzyme (#21)
% and metabolite (#7) as precursor

J{Jcount} = ind_assignV( MS1, [21 7], 1 );

% (4b) C1-C4 (type MS2 = 2) 
% These are catabolic fluxes that recharges ADP to generate ATP. 
% through consuming metabolite (#7-10). It also depends on enzyme (#22)
% and the cencentration of ADP (#12)

J{Jcount} = ind_assignV( MS2, [22 7 12], 2 );
J{Jcount} = ind_assignV( MS2, [22 8 12], 2 );
J{Jcount} = ind_assignV( MS2, [22 9 12], 2 );
J{Jcount} = ind_assignV( MS2, [22 10 12], 2 );

% (5) Define polymerization reaction (type Poly = 4)
% There are 12 polymers and hence 12 synthesis fluxes (Jp1~Jp3, Jq1~Jq8, JR). 
% These polymer synthesis fluxes are proportional through proteome partition ratio,  
% and hence is desbribed by a single synthesis flux.

% The nonlinear synthesis flux is complex and nonlinear, depends on 
% ribosome (#23), ATP (#11), polymer template such as RNA (#24),
% and small metabolite (#7-10) for polymerization.

temp = [23 11 24 7 8 9 10];
J{Jcount} = ind_assignV( Mpoly, temp, 4 );


% (6a) Degradation flux for polymers (#13-#24), (type CD = 5) 
% DP1~DP3, DQ1~DQ8, DR, total 12 fluxes

J{Jcount} = ind_assignV( CD, 13, 5 );
J{Jcount} = ind_assignV( CD, 14, 5 );
J{Jcount} = ind_assignV( CD, 15, 5 );

J{Jcount} = ind_assignV( CD, 16, 5 );
J{Jcount} = ind_assignV( CD, 17, 5 );
J{Jcount} = ind_assignV( CD, 18, 5 );
J{Jcount} = ind_assignV( CD, 19, 5 );
J{Jcount} = ind_assignV( CD, 20, 5 );

J{Jcount} = ind_assignV( CD, 21, 5 );
J{Jcount} = ind_assignV( CD, 22, 5 );

J{Jcount} = ind_assignV( CD, 23, 5 );
J{Jcount} = ind_assignV( CD, 24, 5 );

% (6b) Degradation flux for metabolites (#1-#12), (type CD = 5) 
% Jm1~Jm10, JATP, JADP, total 12 fluxes

J{Jcount} = ind_assignV( CD, 1, 5 );
J{Jcount} = ind_assignV( CD, 2, 5 );
J{Jcount} = ind_assignV( CD, 3, 5 );
J{Jcount} = ind_assignV( CD, 4, 5 );
J{Jcount} = ind_assignV( CD, 5, 5 );
J{Jcount} = ind_assignV( CD, 6, 5 );

J{Jcount} = ind_assignV( CD, 7, 5 );
J{Jcount} = ind_assignV( CD, 8, 5 );
J{Jcount} = ind_assignV( CD, 9, 5 );
J{Jcount} = ind_assignV( CD, 10, 5 );
J{Jcount} = ind_assignV( CD, 11, 5 );
J{Jcount} = ind_assignV( CD, 12, 5 );

%%%

output = param;

% End of the script

