function [ outstring ] = Spectral_determination_script()
%This is a test function for loading scripts into ForestFire. The code may
%look complex because I was testing a complex 'inputdlg' function, but all
%you need to do is just have your custom function output a string that will
%be executed by the cluster. 
%For example: outstring='disp(''''hello world'''');' would work just fine.
%   Philip Andresen, 8-22-2012 (august)

prompt={'Save data?';'Display plot?';'Minimum spectrum limit:';...
    'Spectrum increment:';'Maximum spectrum limit:';...
    'Pixel size (microns)';'Left image boundary (pixels)';...
    'Right image boundary (pixels)';'Bottom image boundary (pixels)';...
    'top image broundary (pixels)';...
    'Calibration file';'Pass previous inputs?'};
name='Spectral Determination';
formats(1,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
formats(2,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
formats(3,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(4,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(5,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(6,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(7,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(8,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(9,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(10,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(11,1)=struct('type','edit','format','file','limits',[0 1],'style','edit');
formats(12,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
defaultanswer={1,1,500,10,820,0.125,1,512,1,512,cd,1};

[inputs,canceled]=inputsdlg(prompt,name,formats,defaultanswer);
if canceled==1; outstring=[]; return; end;
file='original_files';
if inputs{12}==1; file='output_files'; end;
outstring={['Spectral_determination_v04(' ...
    num2str(inputs{1}) ',' num2str(inputs{2}) ',' num2str(inputs{3}) ','...
    num2str(inputs{4}) ',' num2str(inputs{5}) ',' num2str(inputs{6}) ','...
    num2str(inputs{7}) ',' num2str(inputs{8}) ',' num2str(inputs{9}) ','...
    num2str(inputs{10}) ',' '''' inputs{11} '''' ',' file ');']};

end

