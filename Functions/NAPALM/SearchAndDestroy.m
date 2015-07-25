function [ pixel_rcfi_sorted ] = SearchAndDestroy(IMAGE_merged,min_thresh,elim_width)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%elim_width=rbox+rbox2
High_pixel_mask=IMAGE_merged>=min_thresh;
IMAGE_merged=double(IMAGE_merged).*High_pixel_mask;
clear High_pixel_mask
High_pixels_found=uint8(zeros(size(IMAGE_merged)));

convolution_function=ones([2*elim_width+1 2*elim_width+1]);
destruction_matrix=zeros(size(IMAGE_merged));
runs=0;
frames_to_conquer=0;
while max(max(max(IMAGE_merged)))>0;
%     if labindex==1; waitbar((length(frames_to_conquer))/size(IMAGE_merged,3),Loadbar); end;
    requirement=max(max(IMAGE_merged));
    reqmatrix=repmat(requirement,[size(IMAGE_merged,1) size(IMAGE_merged,2) 1]); %the array gives the dimensions of width depth and height
    new_highs=uint8(IMAGE_merged.*(IMAGE_merged==reqmatrix));
    clear reqmatrix requirement
    High_pixels_found=High_pixels_found+(new_highs);
    %destroy the surrounding areas
    frames_to_conquer=find(max(max(IMAGE_merged))>0);
    for i=1:length(frames_to_conquer)
        ind=frames_to_conquer(i);
        destruction_matrix(:,:,ind)=conv2(double(new_highs(:,:,ind)>0),double(convolution_function),'same');
%         if mod(i,20)==0; 
%             waitbar((length(frames_to_conquer))/size(IMAGE_merged,3),Loadbar,...
%                 ['Eliminating adjacent pixels on frames (' num2str(ind) ')']);
%         end;
    end;
    IMAGE_merged(destruction_matrix>0)=0;
    runs=runs+1;
    if runs>5000; 
        disp('Over Run!')
        break; 
    end;
end;

[row_0 col_0]=find(High_pixels_found>0);
Sz=size(High_pixels_found);
pix_index=sub2ind(Sz,row_0,col_0);
[pix_row pix_col pix_frame]=ind2sub(Sz,pix_index);
pix_intensity=double(High_pixels_found(pix_index));
pixel_rcfi=[pix_row,pix_col,pix_frame,pix_intensity]; %row col fram intens
pixel_rcfi_sorted=sortrows(pixel_rcfi,[-4 3]);

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