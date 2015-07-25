function yfit=gaussian_merge3(beta,xdummy)
  global xpix ypix wbox;
 
z=beta(3)*exp(-2*((xpix-beta(1)).^2+(ypix-beta(2)).^2)/(beta(4))^2)+beta(5);
%Z is the intital guess of the gaussan that we will fit!
for i=1:wbox
    for j=1:wbox
      k=(i-1)*wbox+j; %elements [1:j 1:j 1:j 1:j ...] appended along rows
      %Turns the indecies of an x-y grid into a row by row by row list
      yfit(k)=z(i,j);
    end
end

