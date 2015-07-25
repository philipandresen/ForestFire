function yfit=gaussian_merge4(beta,xdummy)
%disp('1')
%rbox=round(beta(6));
wbox=sqrt(size(xdummy,2));
rbox=(wbox-1)/2;
[xpix,ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
%disp('2')
%z=beta(3)*exp(-2*((xpix-beta(1)).*(xpix-beta(1))+(ypix-beta(2)).*(ypix-beta(2)))/(beta(4))^2)+beta(5);
z=beta(3)*exp(-2*((xpix-beta(1)).^2+(ypix-beta(2)).^2)/(beta(4))^2)+beta(5);
%disp('3')
for i=1:wbox
    %disp('4')
    for j=1:wbox
        %disp('5')
      k=(i-1)*wbox+j;
      yfit(k)=double(z(i,j));
    end
end
%disp('6')
