
% Perform inversion and shift image intensity scale

function [output] = Img_invert_shift(input, invert, shift)

input = A1;

invert= 1;
shift = 1000;

output = input;

if (invert == 1)    
    output = (-1)* output;
end

output = output + shift;

figure;
subplot(121);  imagesc(input);
subplot(122);  imagesc(output)

