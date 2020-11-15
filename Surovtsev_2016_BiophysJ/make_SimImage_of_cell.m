function [cell_xy,SimIm3] = make_SimImage_of_cell(n_cells,file_name)
%--------------------------------------------------------------------------
% Simulates fluorescent image of n rectangular cells with m fluorescent molecules.
% Fluorescent molecules are distributed according to chosen distribution
% NOTE: tested only for some distributions
% Meant for test and illustration pruproses 
%--------------------------------------------------------------------------
%
%**********INPUT********:
% n_cells - number of cells 
% file_name - ROI in Matlab sense, nx2 matrix with x znd y coordinates of n vertices of closed ROI
%
% important intrenal variables:
% dx - pixel size
% l0, w0 - cell dimensions
% n_spots - number of lfuorescent molecules
% rand_choice - choice of random distribution rand_XY function
% image_dim - dimensions of simulated image
%
%*********OUTPUT********:
% cell_xy - (x,y,theta)-values for each cell (x,y)-position and theta-oreination angle in image coordinates.
% SimIm3 - simulated image
%
%*********NOTE********:
% Used in Surovtsev, Lim and Jacobs-Wagner. 2016 Biophys.J. v.10 p.2790 
%
%@author:  Ivan Surovtsev
%@date:    2017/10/07
%
%


% cells dimensions in px
dx=0.0642; % px size, um
w0=(0.65-0.1)/dx; % cell length
l0=(2.3-0.1)/dx;  % cell width
l00=l0; w00=w0; 
x_B=0.1*l00;

rand_choice=202; % choice of random distribution

image_dim=[500,500]; %dimensions of simulated image
n_spots0=100; % nuber of spots(molecules)
 %n_cells=64;

% values corresponding to background and single moleucle intensity, based on Hoong's measurements for eLife paper 
bckgrnd=5000;
 bckgrnd_noise=0.014; % relative to backgrnd
signal_height=500;
 signal_noise=0.044;
 signal_width=2.5;

% coefficeint required for distribution according to theory in the BJ paper  
if rand_choice==202,       
  a =  1/(1+0.03*2.3^2/3/0.005);      % based on expressiona_PC/a_tot =  1/(1+k_hyd*l^2/3/D)
  disp(['Fraction of PC-bound ParA=',num2str(a,3)])
end  

% intializing matrices
XYx=repmat((1:image_dim(1)),image_dim(2),1);
 XYy=repmat((1:image_dim(2))',1,image_dim(1));
SimIm=uint16(round(bckgrnd*(1+bckgrnd_noise*randn(image_dim(2),image_dim(1)))));
SimIm3=SimIm;

% simulating cells one by one
for ccell=1:n_cells
  n_spots=round(n_spots0+0.25*n_spots0*randn);    
  x0=(image_dim(1)-2.2*l0)*rand+1.1*l0;
   y0=(image_dim(2)-2.2*l0)*rand+1.1*l0;
  orient=2*pi*rand; % random orientation of the cell in the image
  cell_xy(ccell,:)=[x0,y0,orient,n_spots]; % to store cell position-orientation

  XY0=get_rand_XY(n_spots,rand_choice)'; % get random spots positions
  % rotate cell
  XY=[XY0(:,1)*cos(orient)-XY0(:,2)*sin(orient),XY0(:,1)*sin(orient)+XY0(:,2)*cos(orient)]+repmat([x0,y0],n_spots,1);
   %XY=XY0+repmat([x0,y0],n_spots,1);
  disp(['cell #', num2str(ccell),' coord: ', num2str(x0,3),' ',num2str(y0,3),' totA=', num2str(n_spots)])

% simulate images by sequentially adding image of each cell to the background  
for spot=1:n_spots
  SimIm2=uint16(round(signal_height*exp( -( (XYx-XY(spot,1)).^2+(XYy-XY(spot,2)).^2 )/(signal_width)^2).*(1+abs(signal_noise*randn(image_dim(2),image_dim(1))))));
  SimIm3=SimIm3+SimIm2;
  TotSignal(spot)=sum(SimIm2(:));
end

end


figure; 
 imshow(SimIm3,[min(SimIm3(:)),max(SimIm3(:))],'Border','tight')
 hold on
 %plot(cell_xy(:,1),cell_xy(:,2),'xr')
%figure
 %imshow(SimIm,[min(SimIm(:)),max(SimIm(:))])
 
 
 if nargin>1
   if ischar(file_name)
     % imwrite(SimIm3,file_name,'ColorSpace','cielab');
     imwrite(SimIm3,file_name);
   else
     disp('Can not save the image: 2-nd input "file_name" is not of right type')
   end
 end

function rand_XY=get_rand_XY(nn,choice)
% returns n random XY pairs ([x1,x2,x3...; y1,y2,y3...]) using dPdx according to choice
  uni_rand=rand(2,nn); % uniform random numbers
  switch choice
    case 101 % X: linear gradient from 0 to 1 between current ParB location and new pole
             % Y: uniform distribution
      rand_XY(1,:)=sqrt(uni_rand(1,:))*(l00-x_B)+x_B;
      rand_XY(2,:)=w0*uni_rand(2,:);
    case 201 % X: diffusion/sink distribution from 0 to 1 between current ParB location and new pole
             % Y: uniform distribution
      n_if=1000; % # of points for inverse function approximation
       XXe=0:1/n_if:1; YYe=3*XXe.^2.*(1/2-XXe/6);
      rand_XY(1,:)=(l00-x_B)*interp1(YYe,XXe,uni_rand(1,:),'spline')+x_B; 
       rand_XY(2,:)=w0*uni_rand(2,:);
    case 202 % X: diffusion/sink distribution from 0 to 1 between current ParB location and new pole
             % and accumulation at 0 
             % Y: uniform distribution
      %w0=0.6;
      uni_rand_2=rand(1,nn);
      n_if=1000; % # of points for inverse function approximation
       XXe=0:1/n_if:1; YYe=3*XXe.^2.*(1/2-XXe/6);
       rand_grad=interp1(YYe,XXe,uni_rand(1,:),'spline');
      rand_XY(1,:)=(l00-x_B)*(heaviside(a-uni_rand_2).*zeros(size(rand_grad))+heaviside(uni_rand_2-a).*rand_grad)+x_B; 
       rand_XY(2,:)=w0*(0.5*heaviside(a-uni_rand_2).*ones(size(rand_grad))+heaviside(uni_rand_2-a).*uni_rand(2,:));   
    otherwise
      disp('unknown choice... sorry')  
      return  
  end

end




end

