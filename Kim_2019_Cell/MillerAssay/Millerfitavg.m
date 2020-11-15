function fit = Millerfitavg(X)
%{
-About-
This function averages least-square fit results for Miller data analysis.
This was used to determine time when the first LacZ proteins appear.

-Inputs-
X: (a,b,c) of least-square fit. y = a+bx and y = c 

-Outputs-
fit: fit average

-Example-
fit = Millerfitavg(X);

-Supplementary-

-Keywords-
LacZ, Miller assay

-Dependencies-
millerWorkflow

-References-

-Author-
Sangjin Kim, 2018 April 25
%}

%Define x axis. I used t = 0 to 4 min
x = 0:0.1:4; x = x';
y = x;

%calculate y = a+bx-c for range of x;
%add y = 0 if y falls below 0 (before LacZ appearance)
for i = 1: size(X,1)
    y(:,i) = X(i,2)*x+X(i,1)-X(i,3);
    y(find(y(:,i)<0),i) = 0;
end;

%Average y from different days.
for i = 1:size(y,1)
    yavg(i,1) = mean(y(i,1:end));
    ystd(i,1) = std(y(i,1:end));
end;

figure, 
plot(x,yavg,'o-');

fit = [x, yavg, ystd];
