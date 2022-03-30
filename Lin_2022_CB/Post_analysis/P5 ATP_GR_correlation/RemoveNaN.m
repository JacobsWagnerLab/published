
function [Output] = RemoveNaN(Input, type)

% Remove rows or columns having NaN entry.
%Input = BacDataTotal;
%type = 2;   %1 = row, 2 = column


% (1) flip the matrix if type==2

if (type == 1)  
    
    Mat = Input;    

elseif (type == 2)    

    Mat = Input';    

end


% (2) Looking for rows that contains the flag

rn = size(Mat,1);
temp = isnan(Mat);
 
MarkVector = sum(temp, 2);  
Mat2 = [];

for j = 1:rn
    
    if ( MarkVector(j) == 0 )
        
        Mat2 = [Mat2; Mat(j,:)];
        
    end
    
end


%  (3) flip the matrix if type==2

if (type == 1)  
    
    Output = Mat2;    

elseif (type == 2)    

    Output = Mat2';    

end

