function [outstring]=string_chunkify(instring,maxwidth)
%This little operation takes a string of one line and separates it into 
%many lines of length specified by the limiting factor "maxwidth". 
%words that are too long to fit within maxwidth appear on their own line.
%Due to the nature of 'while' loops it would be advisable to call this 
%function only once per string and save it to a new variable. New lines can
%be dictated with the '//' string and a tab can be dictated with the '/='
%string.
	%Written by Philip Andresen on 7/21/2011.

if iscell(instring);
    for i=1:length(instring)
        if i==1; 
            newinstring=instring{i}; 
        else
            newinstring=[newinstring ' ' instring{i}];
        end;
    end
    instring=newinstring;
end;
        
line=1;
outstring{line}='';
currentline=[];
justbroken=0;
while ~isempty(instring);
	[currentword instring]=strtok(instring);
    switch currentword
        case '//' %If the word is a line break
            line=line+1;
            currentline=[];
            currentword=[];
            justbroken=1;
        case '/=' %if the word is a tab
            currentword='    ';
            justbroken=1;
    end;
    if justbroken==0; currentline=[currentline ' ' currentword]; %#ok<*AGROW>
    else currentline=currentword; justbroken=0; end;
    if length(currentline)<maxwidth;
        outstring{line}=currentline;
    else
        line=line+1;
        currentline=currentword;
        outstring{line}=currentline;
    end;
end
outstring{line}=currentline;
outstring=stack(outstring);
%Stack is an outside function that takes a cell array and makes
%sure it is stacked vertically with whitespace removed.