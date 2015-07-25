function [directory basename extension] = filename_split(full_filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
directory_regexp='\\';
YO=regexp(full_filename,directory_regexp,'ignorecase');
SUP=regexp(full_filename,'\.','ignorecase');
directory=full_filename(1:max(YO));
basename=full_filename(max(YO)+1:SUP-1);
extension=full_filename(SUP:end);
end

