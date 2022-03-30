
function [output] = shuffle_data(input)

L = length(input);
shuffled_index = randperm(L);

output = 0*input;

for j = 1:L        
    output(j) = input( shuffled_index(j) );
end

end
