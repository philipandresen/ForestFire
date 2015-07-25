function image_out=mask_generator(file_in)
load(file_in)
scale_factor=q/0.005;
xf=round(xf_all*scale_factor);
yf=round(yf_all*scale_factor);
[b m n]=unique([xf yf],'rows');
%xw=max(b(:,1));
%yw=max(b(:,2));
image_01=zeros(ceil(xw*scale_factor),ceil(yw*scale_factor));
for i=1:size(b,1)
    image_01(b(i,1),b(i,2))=1;
end;
se=strel('square',ceil(15*scale_factor));
image_02=imopen(image_01,se);
image_out=imclose(image_02,se);