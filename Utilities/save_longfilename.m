function [] = save_longfilename(longfilename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[dir name ext]=filename_split(longfilename);
save([name ext])
copyfile([name ext],dir)

end

