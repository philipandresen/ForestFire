function [ cellarray ] = pullfromstack( cellarray, index )
%PULLFROMSTACK removes an element from an array and restacks the array.
%   Syntax: pullfromstack(array name,index to remove)
%   Written by Philip Andresen
stacksize=size(cellarray);
if stacksize(1,1)<index; index=stacksize(1,1); end;
if isempty(cellarray); return; end;
if stacksize(1,1)==1; cellarray=cellarray'; end
elementnumb=stacksize(1,1);
if index>elementnumb; index=elementnumb; end;
if index<1; index=1; end;
cellarray(index)=[];
cellarray=stack(cellarray);
if size(cellarray)==0; cellarray={''}; end;
%disp(size(cellarray))
end

