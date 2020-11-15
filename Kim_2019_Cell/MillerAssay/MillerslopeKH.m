function mslope = MillerslopeKH(rawData,method)
%{
-About-
This function calculates slope of OD420 increase in Miller assay.

-Inputs-
rawData:  raw data from microplate reader. First column should be time of
reading and columns from second to end should be OD420 readings in each
well.
method: method to select region for fitting. 1 = manual clicking (different
for each well), 2 = define time region (same for all wells)

-Outputs-
mslope:     slope of line fit for each ONPG OD420 increase

-Example-
mslope = MillerslopeKH(rawdata,2);

-Supplementary-

-Keywords-
LacZ, Miller assay

-Dependencies-
millerWorkflow

-References-
Kuhlman, et al. PNAS (2007) vol. 104 pp.6043-8.
https://www.ncbi.nlm.nih.gov/pubmed/17376875

-Author-
Sangjin Kim, 2013 August 02
%}

timeseries = rawData(:,1);
X = rawData(:,2:end);
[timeall,wellnum] = size(X);

% Usual timeseries;
timeseries = [9;19;29;39;49;59;69;79;89;99;109;119;129;149;169;189;209;229;249;269;289;309;329;349;369;389;409;429;449;469;489;509;529;549;569;589;609;629;649;669;689;709;729;749;769;789;809;829;849;869;];

if method == 1
    for i = 1:wellnum 
        display(i)

        % Do manual clinking to select t1 and t2, such that region between t1
        % and t2 is used for line fit
        figure(1), 
        title(['well = ', num2str(i), 'click t1 and t2 for linear regression']);
        plot(timeseries(1:timeall),X(1:timeall,i));hold on;
        pause;
        clk(1:2,1:2) = ginput(2);
        a1 = find(abs(clk(1,1)-timeseries)==min(abs(clk(1,1)-timeseries)));
        a2 = find(abs(clk(2,1)-timeseries)==min(abs(clk(2,1)-timeseries)));
        plot(clk(1,1),clk(1,2),'or',clk(2,1),clk(2,2),'og');
        fit = polyfit(timeseries(a1:a2),X(a1:a2,i),1);
        mslope(i,1) = fit(1);
        hold off;

        % this figure is to check the fit is good
        figure(2),
        plot(timeseries(1:timeall),X(1:timeall,i));hold on;
        plot(timeseries, timeseries.*fit(1)+fit(2),'g-');
        hold off;
    end;
        
elseif method == 2
    figure; cla;
    cc = jet(wellnum);
    for ii = 1:wellnum
        plot(timeseries(1:timeall),X(1:timeall,ii),'Color',cc(ii,:)); hold on;
    end;
    pause;
    clk(1:2,1:2) = ginput(2);
    a1 = find(abs(clk(1,1)-timeseries)==min(abs(clk(1,1)-timeseries)));
    a2 = find(abs(clk(2,1)-timeseries)==min(abs(clk(2,1)-timeseries)));
    plot(clk(1,1),clk(1,2),'or',clk(2,1),clk(2,2),'og');
    hold off;
    
    for i = 1:wellnum
        fit = polyfit(timeseries(a1:a2),X(a1:a2,i),1);
        mslope(i,1) = fit(1);
    end;
    
else
    display('type in correct method as an input');
    return;
end;

