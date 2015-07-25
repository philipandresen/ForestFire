function [ cellarray ] = insert2stack( cellarray, newcontent, index )
%insert2stack will push an element into the stack and shift all elements
%afterwards up one index
%   This will NOT restack the cell array after the insert to prevent
%   accidental whitespace deletion, and the new content is inserted in the 
%   space AFTER the input index.
%   Syntax is: insert2stack(Cell array, new content, index) 
%   Written by Philip Andresen
stacksize=size(cellarray);
newcontentsize=size(newcontent,1);
if stacksize(1,1)==1; cellarray=cellarray'; end
elementnumb=stacksize(1,1);
cellarray{elementnumb+newcontentsize}=[];
if index>elementnumb; index=elementnumb+newcontentsize; end;
if index<1; index=1; end;
for i=index:elementnumb
    cellarray{index+elementnumb-i+newcontentsize}=cellarray{index+elementnumb-i};
end
for i=1:newcontentsize
    cellarray{index+i-1}=newcontent{i};
end;
%cellarray=stack(cellarray);

