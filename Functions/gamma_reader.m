function [X pixel_rcfi_sorted] = gamma_reader(filename)
waitbar_handle=waitbar(0,'Loading, Please wait...');
%J:\Gamma\IMG_1463.mov
%C:\Users\hesslab\Desktop\IMG_1459.mov
obj=VideoReader(filename);
waitbar(0.25,waitbar_handle,'Compiling...')
Video=read(obj);
waitbar(0.99,waitbar_handle,'Preparing for processing...')
hf=figure;
thresh=60;
j=1;
colormap hot
waitbar(0.99,waitbar_handle,'Constructing image stack')
for i=1:obj.NumberOfFrames;
    bwcomposite(:,:,i)=Video(:,:,1,i)+Video(:,:,2,i)+Video(:,:,3,i);
    waitbar(i/obj.NumberOfFrames,waitbar_handle);
end;
clear Video
waitbar(0,waitbar_handle,'Processing...')
for i=1:obj.NumberOfFrames
    floorframe=max(i-20,1);
    ceilframe=min(i+20,obj.NumberOfFrames);
    %bkgn(:,:,i)=mean(bwcomposite(:,:,floorframe:ceilframe),3);
    bkgn=mean(bwcomposite(:,:,floorframe:ceilframe),3);
    meanvalue(i)=mean(mean(bwcomposite(:,:,i)));
    waitbar(i/obj.NumberOfFrames,waitbar_handle);
    bwcomposite(:,:,i)=bwcomposite(:,:,i)-uint8(bkgn);
end;
removed=bwcomposite;
clear bwcomposite bkgn

removed(removed<0)=0;
waitbar(0,waitbar_handle,'Rendering...');
p=1;

j=1;
fh=figure;
ah=axes;
for thresh=116
    %thresh=7*meanvalue;
    specval_sum=0;
    %data=zeros(size(removed));
    disp('Before RCFI')
    pixel_rcfi_sorted=SearchAndDestroy_sequential(removed,thresh,5);
    disp('after_RCFI')
    X(thresh)=size(pixel_rcfi_sorted,1);
    semilogy(X,'Parent',ah);
    disp('After AH')
    waitbar(thresh/250,waitbar_handle)
end;
figure;
    maxval=double(max(max(max(removed))));
    imagesc(maxval-double(max(removed,[],3)))
    for m=1:size(pixel_rcfi_sorted,1);
        hold on
        r=pixel_rcfi_sorted(m,1);
        c=pixel_rcfi_sorted(m,2);
        rectangle('Position',[c-5,r-5,10,10],'EdgeColor',[0 1 0]);
        hold off
    end;
    colormap gray;
    semilogy(X,'Parent',ah)
end

