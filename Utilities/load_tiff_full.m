function [ IMAGE ] = load_tiff_full(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ImageInfo=imfinfo(filename);
num_frames=length(ImageInfo);
waitbarrate=20;
Loadbar=waitbar(0,'Loading image'); %end;
IMAGE=zeros(ImageInfo(1).Height,ImageInfo(1).Width,num_frames); %preallocate memory.
%element=0;
for frame_index=1:num_frames
    %pre-load all selected frames of image.
    IMAGE(:,:,frame_index)=double(imread(filename,frame_index,'Info',ImageInfo));
    %update watibar selectively
    if mod(frame_index,waitbarrate)==0; waitbar(frame_index/num_frames,Loadbar); end;
end

end

