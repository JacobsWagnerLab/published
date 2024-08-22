%% Simple GUI
function pixelshift = pixel_shift_gui(cellList_allfull, dataFilename, brightfield_filename)
close all

global S data h ; %Define this to be global so subfunction can see slider

S.fh = figure;

pixelshiftx = 0.00; pixelshifty = 0.00; 
pixelshift = [pixelshiftx, pixelshifty];

% dataIn = fitsread(brightfield_filename);
dataIn = imread(brightfield_filename);
averageImage = double(mean(dataIn,3));   
    
appData2.data = importdata(dataFilename);
data = appData2.data.data(:,2:3);
clear appData2 dataIn
plotfig(averageImage,cellList_allfull, data, [pixelshiftx pixelshifty]);

hold on
h = plot(data(:,1)+pixelshiftx,data(:,2)+pixelshifty,'b.','MarkerSize',1);


% A slider for varying the parameter.
S.sliderx = uicontrol('Style', 'slider', 'Min',-5,'Max', 5, 'SliderStep',[0.01 0.1],...
'Position', [80 30 200 30], 'Callback', @number_update);

% A slider for varying the parameter.
S.slidery = uicontrol('Style', 'slider', 'Min',-5,'Max', 5, 'SliderStep',[0.01 0.1],...
'Position', [30 60 30 200], 'Callback', @number_update);

S.number = uicontrol('style','edit','units','pix',...
    'position',[350 30 150 40],...
    'string', num2str([pixelshiftx pixelshifty],3),...
    'fontsize',12, 'Callback', @slider_update);

% A button to run the sims.
finish_button = uicontrol('Style', 'pushbutton', 'String', 'Finish',...
'Position', [730 30 100 30], 'Callback', @finish_function);
S.finish = 0;

while S.finish == 0
  pixelshift = str2num(get(S.number,'String'));
  pixelshift = [pixelshift(1) -pixelshift(2)];
  drawnow       
end

close(S.fh);


%% hObject is the button and eventdata is unused.
function number_update(hObject,eventdata)
global data h S;% Slider is a handle to the slider.

% Gets the value of the parameter from the slider.
pixelshiftx = get(S.sliderx,'Value');

% Gets the value of the parameter from the slider.
pixelshifty = get(S.slidery,'Value');

% Puts the value of the parameter on the GUI.
set(S.number,'String', num2str([pixelshiftx pixelshifty],3));

% Plots the Graph.
delete(h); 
h = plot(data(:,1)+pixelshiftx,data(:,2)-pixelshifty,'b.','MarkerSize',1);

function slider_update(hObject,eventdata)
global data h S;% Slider is a handle to the slider.

% Gets the last charactor
last_key_pressed = get(S.fh,'CurrentCharacter');

% If the last_key_pressed is return update the text box
% note that char)13) == Return key
if isequal(last_key_pressed,char(13))
    % This is needed to updated the text stored in this text box
   
   pixelshift = str2num(get(S.number,'String'));
   
   set(S.sliderx,'Value', pixelshift(1))
   set(S.sliderx,'Value', pixelshift(1))
   drawnow
end    

% Plots the Graph.
delete(h); 
h = plot(data(:,1)+pixelshift(1),data(:,2)-pixelshift(2),'b.','MarkerSize',1);

function finish_function(hObject,eventdata)
global data h S finish;% Slider is a handle to the slider.

% Gets the value of the parameter from the slider.
pixelshiftx = get(S.sliderx,'Value');

% Gets the value of the parameter from the slider.
pixelshifty = get(S.slidery,'Value');

% Puts the value of the parameter on the GUI.
S.finish = 1;

% Plots the Graph.
delete(h); 
h = plot(data(:,1)+pixelshiftx,data(:,2)-pixelshifty,'b.','MarkerSize',1);

function plotfig(averageImage,cellList_allfull, data, pixelshift)
imshow(averageImage,[min(min(averageImage)) mean(max(averageImage))], 'InitialMagnification' ,180 );
colormap(gray);
hold on

for ii =  1:length(cellList_allfull)

    
  
    mesh1 = cellList_allfull{ii}.mesh;
    
   
    plgx = [mesh1(1:end-1,1);flipud(mesh1(:,3))];
    plgy = [mesh1(1:end-1,2);flipud(mesh1(:,4))];
    [a,b] = expandpoly(plgx, plgy, 0);  % changed this to 0 from 1
    a(end+1) = a(1);
    b(end+1) = b(1);
    
    xCurve2 = (mesh1(:,1)+mesh1(:,3))/2;
    yCurve2 = (mesh1(:,2)+mesh1(:,4))/2;
        
    % centerline
    plot(xCurve2, yCurve2')
    % outline
    plot(a,b,'color','k');
    axis image
    
end