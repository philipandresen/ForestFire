function [ pixel_rcfi_sorted ] = SearchAndDestroy_sequential(IMAGE_merged,min_thresh,elim_width)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%elim_width=rbox+rbox2
%High_pixel_mask=IMAGE_merged>=min_thresh;
%IMAGE_merged=double(IMAGE_merged).*High_pixel_mask;
IMAGE_merged(IMAGE_merged<min_thresh)=0;
%clear High_pixel_mask
convolution_function=ones([2*elim_width+1 2*elim_width+1]);
imrs=size(IMAGE_merged,1);
imcs=size(IMAGE_merged,2);

tic
pixel_rcfi=[];
numframes=size(IMAGE_merged,3);
hand=waitbar(0,'Finding/eliminating High Pixels');
for i=1:numframes
    currentframe=double(IMAGE_merged(:,:,i));
    while max(max(currentframe))>0
        clear r0 c0 frame0 intens0
        %imagesc(IMAGE_merged(:,:,i));
        %drawnow;
        maxintens=max(max(currentframe));
        newpiximage=currentframe==maxintens;
        [r0 c0]=find(newpiximage>0);
        frame0=ones(size(r0)).*i; %a matrix the same size as r0 but with framenumber
        intens0=ones(size(r0)).*maxintens;
        %new_rcfi=[r0 c0 frame0 intens0]; %RCFI stands for row col frame intens
        pixel_rcfi=[pixel_rcfi;[r0 c0 frame0 intens0]];
        %destruction_image=conv2(single(newpiximage),single(convolution_function),'same');
        newpiximage(max(r0-elim_width,1):min(r0+elim_width,imrs),max(c0-elim_width,1):min(c0+elim_width,imcs))=1;
        %newpiximage=conv2(single(newpiximage),single(convolution_function),'same');
        %destruction_image=newpiximage;
        %logical_destruction_inverse=logical(~destruction_image>0);
        %IMAGE_merged(:,:,i)=IMAGE_merged(:,:,i).*logical_destruction_inverse;
        currentframe=currentframe.*~logical(newpiximage);
    end
    if mod(i,ceil(numframes/400))==0; waitbar(i/numframes,hand); end;
end
%this will prevent an error in sortrows if no elements in pixel_rcfi
toc
close(hand)
if isempty(pixel_rcfi); pixel_rcfi_sorted=[]; return; end;

pixel_rcfi_sorted=sortrows(pixel_rcfi,[-4 3]);
%I've commentend out a lot of code. The uncommented code has been changed 
%to do the same thing as the commented code, but faster (less I/O overhead 
%I guess?) 18.45 seconds down to 4.5 seconds on the benchmark data.


%% Eliminate pixels near the edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this next block of code sets row and column requirements that will
%eliminate pixels that are too close to the edge of the image. First, the
%matrix requirement_matrix contains a repeated list of border coordinates
%and compares those to the rows and columns of our high pixels. Pixels that
%have too short a distance to the edge are deleted.
requirements=[0 0 (size(IMAGE_merged,1)) (size(IMAGE_merged,2))];
requirement_matrix=repmat(requirements,size(pixel_rcfi_sorted,1),1);
Row_col_row_col=[pixel_rcfi_sorted(:,1:2) pixel_rcfi_sorted(:,1:2)];
edge_pixel_distances=Row_col_row_col-requirement_matrix; %compare each row and column to the borders of the image
indecies_to_delete=sum(abs(edge_pixel_distances)<(elim_width+1),2)>0;
pixel_rcfi_sorted(indecies_to_delete,:)=[];%Set each violating element to empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp(['hi ' num2str(length(pixel_rcfi))])
%disp(['hi2 ' num2str(sum(sum(sum(IMAGE_merged>=min_thresh))))])

end