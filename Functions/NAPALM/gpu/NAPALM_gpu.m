%% Fast Localization
% This code is meant to do the same thing as the typical einzel
% localization code, but was written with speed in mind. The code has been
% reorganized into discrete sections, and many operations have been
% modified to exploit the fact that the entire image is loaded into memory
% from the beginning. Recently, the code was fitted with spmd sections or
% 'single program multiple data' sections. This provides parallel
% processing capabilitites.
%% Function name and Batch file GUI help text
function [fname]=NAPALM_gpu(n_start,n_end,bkgn,q,...
    min_thresh,pix_to_pho,rbox,rbox2,rball,corr_file,fname)
%/= func_einzel_2color03 is a program that takes raw data from any image
%file and within each frame attempts to localize molecules based on the
%given parameters. This program is meant for two channel images but can be
%used for more than 2 color imaging depending on how alpha values are
%defined later in the analysis process. // //
%Input variables: //
%/= --------------- //
%/= Start/end frame: //
%/= Select the frames of the input image that should be localized. // //
%/= Background noise: //
%/= Average intensity of background noise in photons per pixel. // //
%/= Pixel size: //
%/= the size of each pixel in the image as measured in micrometers. // //
%/= Minimum Photon Threshold: //
%/= The smallest number of photons detected at a point that is required to
%/= be considered a molecule. This should be in the form of an array with
%/= one value for each file that you process // //
%/= Pixels to photon: //
%/= This is the pixel value that corresponds to one incoming photon. // //
%/= rbox/rbox2: //
%/= rbox is used in a lot of situations. Basically rbox tells us how big of
%a square 'radius' around our molecule we want to capture for analysis. 
%Rbox 2 is used in conjunction to rbox to let the program know when two 
%molecules are too close. // //
%/= rball: //
%/= rball is the radius of the rolling ball usedin the the background 
%subtraction section (slow and complete only)// //
%/=corr_file: //
%/= This is the path and filename of the correlation file to be used in the
%multicolor localization. This correlation file needs to be generated
%before running multicolor localization, and one correlation file may be
%suitable for multiple data sets.
%   This version of the program has been modified to fit with the
%   batch file GUI that is being written for batch file processing. This
%   help file has been modified to work with gethelp.m by Philip Andresen.
tic
clearvars -except n_start n_end bkgn q pix_to_pho rbox rbox2 rball ...
    corr_file fname min_thresh iterations
clc
waitbarrate=5;
%nlabs=4;
%if matlabpool('size')~=0; matlabpool close; end;
%matlabpool('local',nlabs)
%% Load File
% All frames are loaded beforehand to improve speed, but this also has the
% added advantage of freeing up the files for copying, movement, and
% deletion after analysis has moved beyond the loading phase. The 'Imfinfo'
% command is particularly helpful here because when passed to 'imread' it
% allows for continuous reading of the tiff movie.
ImageInfo=imfinfo(fname); %Queue all frames and determine image properties.
if n_end>length(ImageInfo); n_end=length(ImageInfo); end;
if n_start<1; n_start=1; end;
Num_frames=n_end-n_start;
%num_frames=ceil(Num_frames/nlabs);
%framestart=Composite();
framestart=n_start;
frameend=n_end;
num_frames=n_end+1-n_start;
%frameend=Composite();
%for loopindex=1:nlabs
%    framestart{loopindex}=n_start+(loopindex-1)*num_frames;
%    frameend{loopindex}=framestart{loopindex}+num_frames-1;
%    if frameend{loopindex}>n_end;
%        frameend{loopindex}=n_end;
%    end;
%end;

load(corr_file); %For Acoords and Bcoords
if exist('TRANSFORM','var') %If it's a new version correlation file
    Ax1=Acoords(1);
    Ax2=Acoords(3);
    Ay1=Acoords(2);
    Ay2=Acoords(4);
    Bx1=Bcoords(1);
    Bx2=Bcoords(3);
    By1=Bcoords(2);
    By2=Bcoords(4);
    corrtype=1;
else
    load(corr_file,'alpha','beta') %Alpha and beta are bult in matlab functions unless you load them this way
    Ax1=Acoords(1);
    Ax2=Acoords(2);
    Ay1=Acoords(3);
    Ay2=Acoords(4);
    Bx1=Bcoords(1);
    Bx2=Bcoords(2);
    By1=Bcoords(3);
    By2=Bcoords(4);
    corrtype=2;
end;

%spmd
    Loadbar=waitbar(0,'Preloading image'); %end;
    disp(framestart)
    disp(frameend)
    IMAGE=uint16(zeros(ImageInfo(1).Height,ImageInfo(1).Width,num_frames)); %preallocate memory.
    %element=0;
    for frame_index=framestart:frameend
        %element=element+1;
        %pre-load all selected frames of image.
        IMAGE(:,:,frame_index)=uint16(imread(fname,frame_index,'Info',ImageInfo))/uint16(pix_to_pho);
        %update watibar selectively
        if mod(frame_index,waitbarrate)==0; waitbar((frame_index-framestart)/num_frames,Loadbar); end;
    end
    %num_frames=element;
    %% Split channels
    % This code was recently modified to allow for two types of correlation
    % files. The newest correlation program outputs a structure called
    % TRANSFORM that contains a transformation matrix. If that structure is not
    % found we assume that you have selected an old correlation file. A similar
    % if else end statement can be found in the transformation step.
    IMAGE1=IMAGE(Ay1:Ay2,Ax1:Ax2,:); %Channel 1 image stack
    IMAGE2=IMAGE(By1:By2,Bx1:Bx2,:); %channel 2 image stack
%end
clear IMAGE %No need for info outside of region of interest. Save memory!

%% Generate Gaussian
% Before background is removed using the 'slow' method, a gaussian is
% generated that will be used to smooth the noise out of the data. Here is
% where the user can modify the code to select the method of background
% subtraction.
subtraction_type='slow';
if strcmp(subtraction_type,'slow')||strcmp(subtraction_type,'complete')
    FWHM=1.2; %FWHM of gaussian smoothing in pixels
    rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels
    kw=6; %kernal width of smoothing function (Even numbers are best?)
    [Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2);
    kd=sqrt(Xgs.*Xgs+Ygs.*Ygs);
    gs=exp(-2*kd.*kd/(rk*rk));
    gs=gs/sum(sum(gs)); %smoothing function normalized to have area = 1
end;
%spmd
    disp('removing background')
    %% Remove Background
    % There are two methods of background subtraction currently available.
    % Fast, 'mean background' subtraction, and slow 'rolling ball' subtraction.

    % The 'slow' background subtraction type relies on the rolling ball
    % algorithm,  a morphological erosion followed by a morphological
    % dilation using a structuring element that is shaped like a half
    % sphere. Essentially this filter gives you the same result as if
    % you were to track the centre of a large ball as it rolled across
    % a very rough, uneven surface (backgrund+noise) with small holes
    % (data). The holes are not 'detected' by the ball, but the general
    % background levels (large uneven surfaces) are.
    waitbar(0,Loadbar,'Calculating background');
    IMAGE1_noback=zeros(size(IMAGE1,1),size(IMAGE1,2),num_frames); %preallocate memory.
    IMAGE2_noback=zeros(size(IMAGE2,1),size(IMAGE2,2),num_frames); %preallocate memory.
    se = strel('disk',rball,0); %structural element, i.e. rolling ball
    for frame_index=1:num_frames
        %smooth the original and store the frame in 'temp' to use once.
        temp1=uint16(conv2(IMAGE1(:,:,frame_index),gs,'same'));
        temp2=uint16(conv2(IMAGE2(:,:,frame_index),gs,'same'));
        background1 = uint16(imopen(temp1,se));
        background2 = uint16(imopen(temp2,se));
        IMAGE1_noback(:,:,frame_index)=IMAGE1(:,:,frame_index)-background1;
        IMAGE2_noback(:,:,frame_index)=IMAGE2(:,:,frame_index)-background2;
        %update waitbar selectively
        if mod(frame_index,waitbarrate)==0; waitbar((frame_index)/num_frames,Loadbar); end;
    end;

%end
clear IMAGE1 IMAGE2
% if corrtype==1;
%     Temp=Composite();
%     for i=1:nlabs; Temp{i}=TRANSFORM; end;
%     TRANSFORM=Temp;
% else
%     Temp1=Composite();
%     Temp2=Composite();
%     Temp3=Composite();
%     Temp4=Composite();
%     Temp5=Composite();
%     Temp6=Composite();
%     for i=1:nlabs;
%         Temp1{i}=X;  %#ok<*NODEF>
%         Temp2{i}=Y;
%         Temp3{i}=X2;
%         Temp4{i}=Y2;
%         Temp5{i}=alpha;
%         Temp6{i}=beta;
%     end;
%     X=Temp1;
%     Y=Temp2;
%     X2=Temp3;
%     Y2=Temp4;
%     alpha=Temp5;
%     beta=Temp6;
% end;
%spmd
    disp('Transforming Images')
    IMAGE1_noback=IMAGE1_noback.*(IMAGE1_noback>0);
    IMAGE2_noback=IMAGE2_noback.*(IMAGE2_noback>0);
    
    %% Transform channel 2 to match channel 1
    % This code will operate differently depending on the version of
    % correlation file that has been loaded.
    %%
    % The first listed method is for the new correlation file type. A built in
    % matlab function transforms the image using a transformation structure
    % from the correlation file. *NOTE: for this method XData and YData MUST be
    % specified*, otherwise the image will be translated to best fit within the
    % original image size. We need to maintain the same coordinate system
    % throughout analysis, so this added specification is imperative
    %%
    % The second listed method uses a pregenerated set of coordinates to
    % transform the image. In the old correlation system, X and Y were meshgrid
    % generated coordinates that would be transformed into X2 Y2. We can
    % effectively transform the image by using an interpolation function on the
    % image between these two coordinate systems.
    %%
    % After the image is transformed, the program ensures that the total number
    % of photons has remained constant by scaling the new intensity data such
    % that the total intensity in the transformed channel is the same as it was
    % before transformation.
    waitbar(0,Loadbar,'Transforming Channel 2 onto channel 1'); %end;
    IMAGE2_transf_noback=zeros(size(IMAGE2_noback,1),size(IMAGE2_noback,2),num_frames); %preallocate memory.
    for frame_index=1:num_frames
        if corrtype==1;
            A=IMAGE2_noback(:,:,frame_index);
            IMAGE2_transf_noback(:,:,frame_index)=imtransform(A,TRANSFORM,...
                'XData',[1 size(A,2)],'YData',[1 size(A,1)]); %THIS IS SUPER IMPORTANT Read help for imtransform for info.
            integral1=sum(IMAGE2_noback(:));
            IMAGE2_transf_noback(isnan(IMAGE2_transf_noback))=0; %set any 'NaN' to 0
            integral2=sum(IMAGE2_transf_noback(:));
            IMAGE2_transf_noback=IMAGE2_transf_noback.*(integral1/integral2);
            %Make sure you have the same number of photons. (four previous lines)
        else
            IMAGE2_transf_noback(:,:,frame_index)=interp2(X,Y,IMAGE2_noback(:,:,frame_index),X2,Y2)/(alpha*beta);
            IMAGE2_transf_noback(isnan(IMAGE2_transf_noback))=0; %set any 'NaN' to 0
            %Alpha*beta has to do with the jacobian?
        end;
        if mod(frame_index,waitbarrate)==0; waitbar((frame_index)/num_frames,Loadbar); end;
    end;
    
    disp(sum(IMAGE2_transf_noback(:)))
    disp('Done transforming Images')
%end
clear IMAGE2_noback %The untransformed version of IMAGE2 is no longer needed.

%spmd
    disp('Finding potential fluorophores')
    %% Construct a channel merged version and free memory
    IMAGE_merged=IMAGE1_noback+IMAGE2_transf_noback; %Combine images.
    % clearvars -except IMAGE_merged IMAGE2_transf_noback Num_frames min_thresh ...
    %     Loadbar rbox rbox2 q bkgn pix_to_pho ImageInfo fname IMAGE1_noback ...
    %     waitbarrate num_frames
    
    %% Find and eliminate nearby high pixels
    % The SeachAndDestroy algorithm was written to address slow pixel
    % elimination times for files that contained many above-threshold pixels.
    % This algorithm uses 2D convolutions and 3D array operations to eliminate
    % as many adjacent pixels as quickly as possible while making sure that the
    % brightest pixels have priority. Rbox+rbox2 is the width of the
    % elimination square.
    waitbar(0,Loadbar,'Finding and eliminating adjacent high pixels'); %end;
    %pixel_rcfi_sorted=SearchAndDestroy(IMAGE_merged,min_thresh,elim_width);
    pixel_rcfi_sorted=SearchAndDestroy_sequential(IMAGE_merged,min_thresh,rbox+rbox2);
    disp('Fitting Point spread functions')
    %% Fit PSFs
    % This code can be separated into three basic steps.
    %
    % *Grabbing the Data*
    %
    % The program selects a high pixel isolated in the previous section, it then
    % grabs a small image of the pixel's neighborhood which is usually a 7x7
    % square region. This 'grab' is used to guess an x and y using a
    % weighted center of mass calculation.
    %
    % *Data Unfolding*
    %
    % The program creates a 2D
    % gaussian that might fit the grabbed pixels. In order to actually test
    % whether the gaussian fits, it is turned into a 1D function by being
    % sliced up and unfolded into a series of growing and shrinking gaussians.
    % The same is done for the actual grabbed data.
    %
    % *'Nlinfit'*
    %
    % a 1D fitting
    % algorithm called 'nlinfit' tries to fit the unfolded 2D gaussian to the
    % unfolded 2D data. Because both have been unfolded in the same way a 2D
    % fit can be derived from the parameters defined in the 1D 'nlinfit'.
    n_boxes=0;
    box_overlap_factor = 1.5; %if center to center closer, don't include either
    w_mask = round(rbox*box_overlap_factor);
    psf_scale=1.2;
    NA=1.2;
    wvlnth=580/1000; %let the user define this too, please!
    psf_w0 = psf_scale*0.55*wvlnth/NA/1.17;
    psf_std=psf_w0/2; %standard deviation of psf
    psf_w02=(psf_w0/q)*(psf_w0/q);  %square of 1/e^2 radius in pixels
    %global wbox xpix ypix
    wbox=2*rbox+1; %#ok<*NASGU>
    [xpix0,ypix0] = meshgrid(-2*rbox:2*rbox,-2*rbox:2*rbox);  %#ok<NASGU,ASGLU>
    [xpix,ypix] = meshgrid(-rbox:rbox,-rbox:rbox);
    total_molecules=0;
    n_fail_a0=0;
    n_fail_outbox=0;
    %Initialize data arrays
    n_init=length(pixel_rcfi_sorted);
    [xcm_all ycm_all xf_all yf_all a0_all r0_all off_all framenum_all ...
        xf_err_all yf_err_all a0_err_all r0_err_all off_err_all grab_sum_all...
        green_sum red_sum]=deal(zeros(n_init,1));
    boxes_xy=zeros(10000,3);
    %end data inititalization
    waitbar(0,Loadbar,['Fitting Point Spread Functions. (' num2str(size(pixel_rcfi_sorted)) ')']);
    warning('OFF') %#ok<WNOFF>
    iterations=100;
    options=statset('Maxiter',iterations);
    if iterations<100;
        disp(['Be careful! Curve fitting set to only ' num2str(iterations) ' iterations! (default: 100)'])
    end;
    for pixel_index=1:size(pixel_rcfi_sorted)
        %%
        % Here we have the row, column, and frame of the high pixel.
        R=pixel_rcfi_sorted(pixel_index,1); %Current pix row
        C=pixel_rcfi_sorted(pixel_index,2); %Current pix col
        F=pixel_rcfi_sorted(pixel_index,3); %Current pix frame
        %%
        % This is the grabbed pixel neighborhood
        grab=IMAGE_merged((R-rbox):(R+rbox),(C-rbox):(C+rbox),F);
        %%
        % This is the unfolded 2D data and best x-y guesses (by center of mass).
        zmerge = reshape(grab',1,[]);
        if pixel_index==1; xymerge = zeros(1,numel(grab)); end;
        %Center of mass guess for fit (Tested code, works!)
        bestguess_col=sum(sum((grab.*xpix)/sum(grab(:))));
        bestguess_row=sum(sum((grab.*ypix)/sum(grab(:))));
        %%
        % These are the parameters for the unfolded 2D gaussian that we mean to fit!
        Init_guess_param=[bestguess_col,bestguess_row,50,psf_w0/q,min(grab(:)),rbox]; % c,r,a0,r0,offset
        %%
        % This is the actual fitting line.
        [fit_param,resid,~,COVB,~] = nlinfit(xymerge,zmerge,@gaussian_merge4,Init_guess_param,options);
        %%
        confidence_intervals = nlparci(fit_param,resid,'covar',COVB); %calculate error estimates on parameters
        confidence_intervals_err=(confidence_intervals(:,2)-confidence_intervals(:,1))/2;
        y_fit=fit_param(2)+R;
        x_fit=fit_param(1)+C;
        a0_fit=fit_param(3);
        r0_fit=abs(fit_param(4));
        offset=fit_param(5);
        if (x_fit>(C+rbox) || x_fit<(C-rbox) || y_fit<(R-rbox) || y_fit>(R+rbox));
            n_fail_outbox=n_fail_outbox+1;
            continue; %Skip saving if outside of rbox region
        end;
        if (x_fit<(rbox) || x_fit>(size(IMAGE_merged,2)-rbox) || y_fit<(rbox) || y_fit>(size(IMAGE_merged,1)-rbox));
            n_fail_outbox=n_fail_outbox+1;
            continue; %Skip saving if outside of rbox region
        end;
        if a0_fit<0;
            n_fail_a0=n_fail_a0+1;
            continue; %Skip saving if negative amplitude
        end;
        %save
        total_molecules=total_molecules+1;
        xcm_all(total_molecules)=C ;
        ycm_all(total_molecules)=R;
        xf_all(total_molecules)=x_fit;
        yf_all(total_molecules)=y_fit;
        a0_all(total_molecules)=a0_fit;
        r0_all(total_molecules)=r0_fit;
        off_all(total_molecules)=offset;
        framenum_all(total_molecules)=F;
        xf_err_all(total_molecules)=confidence_intervals_err(1);
        yf_err_all(total_molecules)=confidence_intervals_err(2);
        a0_err_all(total_molecules)=confidence_intervals_err(3);
        r0_err_all(total_molecules)=confidence_intervals_err(4);
        off_err_all(total_molecules)=confidence_intervals_err(5);
        grab_sum_all(total_molecules)=-1;
        
        im_green=IMAGE1_noback((round(y_fit)-rbox2):(round(y_fit)+rbox2),(round(x_fit)-rbox2):(round(x_fit)+rbox2),F);
        im_red=IMAGE2_transf_noback((round(y_fit)-rbox2):(round(y_fit)+rbox2),(round(x_fit)-rbox2):(round(x_fit)+rbox2),F);
        green_sum(total_molecules)=sum(im_green(:));
        red_sum(total_molecules)=sum(im_red(:));
        
        n_boxes=n_boxes+1;
        boxes_xy(n_boxes,1)=C;
        boxes_xy(n_boxes,2)=R;
        boxes_xy(n_boxes,3)=1;
        
        if mod(pixel_index,waitbarrate)==0; waitbar(pixel_index/size(pixel_rcfi_sorted,1),Loadbar); end;
    end;
    disp([total_molecules Num_frames n_fail_a0 n_fail_outbox])
%end

% %%%%%%%%%%
% IMAGE_merged2=cat(3,IMAGE_merged{:});
% full_path=cd;
% total_frames=size(IMAGE_merged2,3);
% [output_file, output_path]=uiputfile('*.tif','Choose a location for your file',full_path);
% full_path=[output_path output_file];
% for y=1:1;
%     ops.color   = false;
%     ops.comp    = 'no';
%     ops.message = true;
%     ops.ask     = true;
%     ops.append  = true;
%     if y==1; ops.append=false; end;
%     saveastiff(uint16(IMAGE_merged2(:,:,:)), full_path, ops)
%     %imwrite(stacked_tiff(:,:,y),full_path,'tif','WriteMode','append')
%     waitbar((y-1)/total_frames,waitbar_id,'writing')
% end
% close(waitbar_id)
% %%%%%%%%%



%% Output Completion Stats
% These are simply some aesthetic commands to tell the user how long
% program execution took. Sometimes debug statements will show up here as
% well.
Time_taken=toc;
disp(['Completed in ' num2str(floor(Time_taken/60)) ' minutes and ' num2str(round(mod(Time_taken,60))) ' seconds.'])
%waitbar(1,Loadbar,'Complete','Color',[1 1 1]);
%BarHandle=findobj(Loadbar,'Type','Patch');
%set(BarHandle,'EdgeColor',[0 0 0],'FaceColor',[0 1 0]) %changes the color to green
%spmd
    %% Truncate Arrays and Save
    % All arrays were preallocated to be too long. Now they are all trimmed and
    % saved. This can take a little time, especially since output histograms
    % are generated and saved here.
    xcm_all=xcm_all(1:total_molecules);
    ycm_all=ycm_all(1:total_molecules);
    xf_all=xf_all(1:total_molecules);
    yf_all=yf_all(1:total_molecules);
    a0_all=a0_all(1:total_molecules);
    r0_all=r0_all(1:total_molecules);
    off_all=off_all(1:total_molecules);
    framenum_all=framenum_all(1:total_molecules)+(labindex-1)*num_frames;
    xf_err_all=xf_err_all(1:total_molecules);
    yf_err_all=yf_err_all(1:total_molecules);
    a0_err_all=a0_err_all(1:total_molecules);
    r0_err_all=r0_err_all(1:total_molecules);
    off_err_all=off_err_all(1:total_molecules);
    grab_sum_all=grab_sum_all(1:total_molecules);
    green_sum=green_sum(1:total_molecules);
    red_sum=red_sum(1:total_molecules);
%end
% total_molecules=sum([composite_total_molecules{:}]);
% n_fail_a0=sum([composite_n_fail_a0{:}]);
% n_fail_outbox=sum([composite_n_fail_outbox{:}]);
% xcm_all=vertcat(composite_xcm_all{:});
% ycm_all=vertcat(composite_ycm_all{:});
% xf_all=vertcat(composite_xf_all{:});
% yf_all=vertcat(composite_yf_all{:});
% a0_all=vertcat(composite_a0_all{:});
% r0_all=vertcat(composite_r0_all{:});
% off_all=vertcat(composite_off_all{:});
% framenum_all=vertcat(composite_framenum_all{:});
% xf_err_all=vertcat(composite_xf_err_all{:});
% yf_err_all=vertcat(composite_yf_err_all{:});
% a0_err_all=vertcat(composite_a0_err_all{:});
% r0_err_all=vertcat(composite_r0_err_all{:});
% off_err_all=vertcat(composite_off_err_all{:});
% grab_sum_all=vertcat(composite_grab_sum_all{:});
% green_sum=vertcat(composite_green_sum{:});
% red_sum=vertcat(composite_red_sum{:});
% box_overlap_factor=box_overlap_factor{1};
% X=X{1};
% X2=X2{1};
% Y=Y{1};
% Y2=Y2{1};
% alpha=alpha{1};
% beta=beta{1};
% framestart=vertcat(framestart{:});
% frameend=vertcat(frameend{:});
% iterations=iterations{1};
% n_boxes=sum([n_boxes{:}]);
% psf_scale=psf_scale{1};
% psf_w0=psf_w0{1};
% psf_std=psf_std{1};
% psf_w02=psf_w02{1};
% wvlnth=wvlnth{1};
% NA=NA{1};
% 
% matlabpool close


nrat=red_sum./(red_sum+green_sum);
bins1=0:0.01:1;
h1=hist(nrat,bins1);
ratio=red_sum./green_sum;
bins2=0:0.1:20;
h2=hist(ratio,bins2);
alphagraph=figure;
title(fname);
subplot(2,1,1)
bar(bins2,h2);
xlabel('Red/Green');
ylabel('#');
xlim([0 6]);
subplot(2,1,2)
bar(bins1,h1);
xlabel('Red/(Red+Green)');
ylabel('Frequency');
legend(fname);
xlim([0 1]);
print(alphagraph,'-dbmp',strrep(fname,'.tif','_alphas.bmp'))

npix_all=pi*(r0_all.^2)/2;     % area of molecule in square pixels
N=npix_all.*a0_all; % number of photons for each molecule
lp2=((r0_all*q/2).^2+(q^2)/12)*1./N+8*pi*((r0_all*q/2).^4)*(bkgn^2)/(q^2)*1./(N.*N);
lp=1.3*sqrt(lp2); %loc prec in um (add 30%, Thompson et al, 2002)
locprec=figure;
title(fname);
hist(lp*1000,100);
xlabel('Loc. Prec. (nm)');
ylabel('#');
print(locprec,'-dbmp',strrep(fname,'.tif','_LocPrec.bmp'))

clear IMAGE_merged IMAGE1_noback IMAGE2_transf_noback IMAGE1 IMAGE2_transf ...
    pixel_rcfi pixel_rcfi_sorted
fname=strrep(fname,'.tiff','');
fname=strrep(fname,'.tif','');
fname=[fname ' n=' num2str(total_molecules) '.mat'];

clear COVB Init_guess_param num_frames ...
    R C F I Temp1 Temp2 Temp3 Temp4 Temp5 Temp6 background1 background2 ...
    bestguess_row a0_fit boxes_xy a0_all a0_err_all...
    framenum_all grab_sum_all green_sum...
    n_fail_a0 n_fail_outbox off_all...
    off_err_all r0_all r0_err_all...
    red_sum total_molecules xcm_all...
    xf_all xf_err_all ycm_all...
    yf_all yf_err_all confidence_intervals...
    confidence_intervals_err element fit_param frame_index grab im_green...
    im_red integral1 integral2 n_init offset options pixel_index...
    r0_fit resid se temp1 temp2 w_mask wbox x_fit xpix xpix0 y_fit ypix...
    ypix0 zmerge bestguess_col xymerge M1 M2 Min1 Min2 Background_special1 ...
    Background_special2 extra1 extra2 factor1 factor2
save(fname)
end