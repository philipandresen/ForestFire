function [ outstring ] =select_species_one_species_call3_script()
%This is a test function for loading scripts into ForestFire. The code may
%look complex because I was testing a complex 'inputdlg' function, but all
%you need to do is just have your custom function output a string that will
%be executed by the cluster. 
%For example: outstring='disp(''''hello world'''');' would work just fine.
%   Philip Andresen, 8-22-2012 (august)
%[]=select_species_one_species_call3(um_per_pix,drmax,npmax,min_points_per_cluster,infile)
prompt={'um per pixel';'drmax';'npmax';...
    'minimum points per cluster';'Species? (0=none,1=yellow,2=red)';...
    'Green alpha limits';'Red alpha limits';'Input previous?'};
name='Cluster analysis';
formats(1,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(2,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(3,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(4,1)=struct('type','edit','format','integer','limits',[0 inf],'style','edit');
formats(5,1)=struct('type','edit','format','integer','limits',[0 2],'style','edit');
formats(6,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(7,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(8,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
defaultanswer={0.137;0.03;25000;200;0;'[0 0.55]';'[0.68 1]';1};

[inputs,canceled]=inputsdlg(prompt,name,formats,defaultanswer);
if canceled==1; outstring=[]; return; end;
file='original_files';
if inputs{8}==1; file='output_files'; end;
outstring={['select_species_one_species_call3(' ...
    num2str(inputs{1}) ',' num2str(inputs{2}) ',' num2str(inputs{3}) ','...
    num2str(inputs{4}) ','  num2str(inputs{5}) ',' inputs{6} ',' inputs{7} ',' file ');']};
end