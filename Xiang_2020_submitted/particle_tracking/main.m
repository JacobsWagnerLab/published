%% Parse in all the .tif stacks
fls = dir('*.tif');

%% Parameters array
ths =    [];
minSzs = [];
maxSzs = [];

%% Loop through individual stack
for ii = 1 : length(fls)
%% Construct the SPT object
thisStack = [pwd,'/',fls(ii).name];
fprintf('--------------------\n');
fprintf('It is now %s\n', datestr(now));
fprintf('Working on %s, %d/%d left for analysis... \n', fls(ii).name(1:end-4), length(fls)-ii, length(fls));

date = '0000-00-00';
about = 'Here goes some information about the experiment';
stackPath = thisStack;
pixelLength = 1/9.3; % um/px
particleSize = 1; % um
frameTime = 0.5; % sec

p = spt('date',date,'about',about,'stackPath',stackPath,...
        'particleSize',particleSize,...
        'pixelLength',pixelLength,'frameTime',frameTime);

%% Set the filtering threshold
Ithreshold = ths(ii);
% p.sampleFrame(Ithreshold);

%% Get all particle locations
minSpotSize = minSzs(ii);
maxSpotSize = maxSzs(ii);
eccentricity = 0.5;

tic
p.getLocation('minspotsize',minSpotSize,...
              'maxspotsize',maxSpotSize,...
              'eccentricity', eccentricity, ...
              'intensityThreshold',Ithreshold);                                 
toc

save([fls(ii).name(1:end-4),'.mat']);
clear thisStack p Ithreshold minSpotSize maxSpotSize;

if (ii == 1)
    mailMe('sekainokami@gmail.com','You first analysis is done!','');
end

end
sendTo = 'sekainokami@gmail.com';
subject = 'Analysis done!';
text = ['All analyses are finished at ', datestr(datetime)]; 
mailMe(sendTo,subject,text);

%% Check localization
p.plotLocation();

%% test
tic
p.getLocation('minspotsize',50,...
              'maxspotsize',200,...
              'eccentricity', 0.5, ...
              'intensityThreshold',150,'startframe',1,'endframe',10);                                 
toc