function [ outstring ] = NAPALM_script()
%This is a test function for loading scripts into ForestFire. The code may
%look complex because I was testing a complex 'inputdlg' function, but all
%you need to do is just have your custom function output a string that will
%be executed by the cluster. 
%For example: outstring='disp(''''hello world'''');' would work just fine.
%   Philip Andresen, 8-22-2012 (august)

prompt={'First frame';'Last frame';'Background noise (photons)';...
    'Pixel size (microns)';'Minimum threshold (photons)';...
    'Pixel to photon ratio';'rbox';'rbox2';'rball';'correlation file';...
    'Input the previous output?'};
name='NAPALM initialization';
formats(1,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(2,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(3,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(4,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(5,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(6,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(7,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(8,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(9,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(10,1)=struct('type','edit','format','file','limits',[0 1],'style','edit');
formats(11,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
defaultanswer={1;10000;4.85;0.136;'[20]';40;3;2;6;cd;1};

[inputs,canceled]=inputsdlg(prompt,name,formats,defaultanswer);
if canceled==1; outstring=[]; return; end;
file='original_files';
if inputs{11}==1; file='output_files'; end;
outstring={['threshold= ' inputs{5} ';'];...
    ['output_files = NAPALM(' ...
    num2str(inputs{1}) ',' num2str(inputs{2}) ',' num2str(inputs{3}) ','...
    num2str(inputs{4}) ',' 'threshold(currentfile)' ',' num2str(inputs{6}) ','...
    num2str(inputs{7}) ',' num2str(inputs{8}) ',' num2str(inputs{9}) ','...
    '''' inputs{10} '''' ',' file ');']};

end

