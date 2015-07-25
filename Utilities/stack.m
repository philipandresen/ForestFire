function [ cellarray ] = stack( cellarray )
%Stack creates a stack that can have elements inserted into it and elements
%removed from it while keeping the cell array tightly stacked
%   Calling stack on a cell/numeric array will cause the elements to be stacked
%   vertically with all whitespace removed. It will NOT remove the last
%   item in a list if it is whitespace. Calling stack on a character type
%   variable or a non-array will cause an error to prevent verticalization 
%   of strings which is often undesireable and is not the purpose of this
%   program.
%   Written by Philip Andresen
if iscell(cellarray)==0; 
    if isnumeric(cellarray); 
        cellarray=num2cell(cellarray);
        disp('WARNING! Numeric array converted to Cell array!')
    else
        disp('Input Argument must be Cell array, X{a b c...}'); 
        return; 
    end
end;

stacksize=size(cellarray);
if stacksize(1,1)==1; cellarray=cellarray'; end
if iscellstr(cellarray)==1; cellarray=deblank(cellarray); end
stacksize=size(cellarray);
elementnumb=stacksize(1,1);
for i=1:elementnumb-1
    switch ischar(cellarray{elementnumb-i})
        case 1
            if strcmp(cellarray(elementnumb-i),'')==1; %If blank, shift all later elements down
                cellarray(elementnumb-i)=[];
            end
        case 0
            if cellarray{elementnumb-i}==0;
                cellarray(elementnumb-i)=[];
            end
    end
end
