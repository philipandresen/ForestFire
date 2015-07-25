function [ outstring ] = Einzel_script_prism()
%This is a test function for loading scripts into ForestFire. The code may
%look complex because I was testing a complex 'inputdlg' function, but all
%you need to do is just have your custom function output a string that will
%be executed by the cluster. 
%For example: outstring='disp(''''hello world'''');' would work just fine.
%   Philip Andresen, 8-22-2012 (august)

prompt={'Pixel size (microns)';'Minimum photon threshold';'Pixel to photon conversion';...
    'PSF region size (rbox)';'Background subtraction radius (rball)';...
    'Max particles per frame';'Left image boundary (pixels)';...
    'Right image boundary (pixels)';'Bottom image boundary (pixels)';...
    'top image broundary (pixels)';...
    'Calibration file'};
name='Einzel initialization';
formats(1,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit'); %
formats(2,1)=struct('type','edit','format','integer','limits',[0 inf],'style','edit');%
formats(3,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit'); %
formats(4,1)=struct('type','edit','format','integer','limits',[0 inf],'style','edit'); %
formats(5,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit'); %
formats(6,1)=struct('type','edit','format','integer','limits',[0 inf],'style','edit'); %
formats(7,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit'); %
formats(8,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit'); %
formats(9,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit'); %
formats(10,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit'); %
formats(11,1)=struct('type','edit','format','file','limits',[0 1],'style','edit'); %
defaultanswer={0.12;10;20;3;12;15;1;512;1;512;cd};

[inputs,canceled]=inputsdlg(prompt,name,formats,defaultanswer);
if canceled==1; outstring=[]; return; end;
file='original_files';
if inputs{11}==1; file='output_files'; end;
outstring={['output_files = func_einzel_prism_MM_v01_multi(' ...
    num2str(inputs{1}) ',' num2str(inputs{2}) ',' num2str(inputs{3}) ','...
    num2str(inputs{4}) ',' num2str(inputs{5}) ',' num2str(inputs{6}) ','...
    num2str(inputs{7}) ',' num2str(inputs{8}) ',' num2str(inputs{9}) ','...
    num2str(inputs{10}) ',' file ',' '''' inputs{11} '''' ');']};

end

