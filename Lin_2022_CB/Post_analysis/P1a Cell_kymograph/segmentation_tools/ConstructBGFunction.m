function [bg_function, bg_matrix] = ConstructBGFunction(image,seg_mat,sampling_interval,plot)
%{
-About-
This function reconstrucs the background of an image (fluorescente or
phase) with a two-dimensional polynome of second degree. The resulting
function can be used to get a good estimate of the background at a certain
position (pixel in image). As the Background often is not uniform this
allowes accurate local background subtraction over the whole image.
Note: This mostly applies for images of agarose-pads

-Inputs-
image:              An n by m matrix representing an image.
seg_mat:            An n by m matrix where non-zero-values represent the
                    pixel mask of objects (eg. non-background) present.
sampling_interval:  (optional) An Integrer which defines the interval
                    between two saples taken for the approximation.
                    Default = size(image,1)/64 (eg. 32 for 2048x2048)
plot:               (optional) A logical which when true allowes to plot
                    the samples and the reconstruction for visualisation.
                    Default = false

-varargin-
none

-Outputs-
bg_function:     A function F(x Col,y Row) which calculates the background
                 value for the desired point.

-Example-
ConstructBGFunction(image_channel_3,cell_mask,16,true)
   
-Supplementary-
Supplemental file location

-Keywords-
Background, Construction, Reconstruction, subtraction, agarose-pad, pad,
agarose,

-Dependencies-

-References-
None

-Author-
Nikolaus Huwiler, 2017 August 25
%}



% Ensure robustness of input arguments
if nargin < 4
    plot = false;
end
if nargin < 3
    sampling_interval = size(image,1)/64;
end
if nargin < 2
    error('Not enough input arguments');
end


% allocate the matrix which is chose to select the samples based on the
% chosen sample interval
bg_samples = false(size(seg_mat));

for i = 1:size(bg_samples,1)
    for j = 1:size(bg_samples,2)
        if not(logical(rem(i,sampling_interval))) && not(logical(rem(j,sampling_interval)))
            bg_samples(i,j) = true;
        end
    end
end

% Ensure that no pixel is taken into account which is not part of the
% background
bg_samples = logical(bg_samples .* not(logical(seg_mat)));


% Now that the position of the samples have been selected we asociate the
% values of the image to them
bg_samples_for_approx = zeros(size(bg_samples));
bg_samples_for_approx(bg_samples) = image(bg_samples);

% Make the polynomial approximation with the samples chosen and build the
% output fuction
[y,x,z] = find(bg_samples_for_approx);
f = fit( [x, y], z, 'poly22' );
p = coeffvalues(f); %a + b*x + c*y + d*x^2 + e*x*y + f*y^2
bg_function = @(c,r) p(1) + p(2).*c + p(3).*r + p(4).*c.^2 + p(5).*c.*r + p(6).*r.^2;

bg_matrix = 0*image;

for r = 1:size(image,1)
    for c = 1:size(image,2)
        bg_matrix(r,c) = bg_function(r,c);
    end
end

% plot the sampes and the function
if plot
    [xx,yy] = meshgrid(1:10:size(bg_samples,1),1:10:size(bg_samples,2));
    figure();
    p1 = surf(xx,yy,bg_function(xx,yy));
    grid off;
    hold on;
    p2 = scatter3(x,y,z,'ro');
    legend([p1(1) p2(1)],{'BG approx','BG samples'});
    xlabel('X');
    ylabel('Y');
end





end