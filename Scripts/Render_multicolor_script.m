function [ outstring ] = Render_multicolor_script()
%This is NOT a test function for loading scripts into ForestFire. The code may
%look complex because I was testing a complex 'inputdlg' function, but all
%you need to do is just have your custom function output a string that will
%be executed by the cluster. 
%For example: outstring='disp(''''hello world'''');' would work just fine.
%   Philip Andresen, 8-22-2012 (august)


prompt={'Save rendered image (y or n)=','Number of Species=',...
             'Species 1 range=','Species 2 range','Species 3 range',...
             'Image scaling factor=','Molecule scaling factor',...
             'Sp.1 RGB disp. color=','Sp. 2 RGB disp. color=',...
             'Sp. 3 RGB disp. color=','Pass previous output to this program?'};
name='Multicolor Render initialization';
formats(1,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(2,1)=struct('type','edit','format','integer','limits',[1 3],'style','edit');
formats(3,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(4,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(5,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(6,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(7,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(8,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(9,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(10,1)=struct('type','edit','format','text','limits',[0 1],'style','edit');
formats(11,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
defaultanswer={'y',1,'[0.3 0.6]','[0.7 0.85]','[0.9 1.0]',14,2.0,'[0 1 0]','[1 0 0]','[0 0.5 1]',1};

[inputs,canceled]=inputsdlg(prompt,name,formats,defaultanswer);
if canceled==1; outstring=[]; return; end;
file='original_files';
if inputs{11}==1; file='output_files'; end;
outstring={['function_render_multicolor(' ...
    '''' inputs{1} '''' ',' num2str(inputs{2}) ',' num2str(inputs{3}) ','...
    num2str(inputs{4}) ',' '[0.85 1]' ',' num2str(inputs{6}) ','...
    num2str(inputs{7}) ',' num2str(inputs{8}) ',' num2str(inputs{9}) ','...
     inputs{10} ',' file ');']};

end

