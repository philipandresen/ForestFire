function [ outstring ] = bleedthrough_correction_script()
%This is a test function for loading scripts into ForestFire. The code may
%look complex because I was testing a complex 'inputdlg' function, but all
%you need to do is just have your custom function output a string that will
%be executed by the cluster. 
%For example: outstring='disp(''''hello world'''');' would work just fine.
%   Philip Andresen, 8-22-2012 (august)
%[eff_d2c eff_c2d d_min d_max c_min c_max]=bleedthru_corr('dendra','Cherry',Dup,Cdn)
prompt={'Dendra file (single species, two channel)';...
    'Cherry file (single species, two channel)';...
    'Upper dendra alpha ratio';'Lower cherry alpha ratio';...
    'Display histogram?'};
name='Bleedthrough Correction';
formats(1,1)=struct('type','edit','format','file','limits',[0 1],'style','edit');
formats(2,1)=struct('type','edit','format','file','limits',[0 1],'style','edit');
formats(3,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(4,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(5,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
defaultanswer={cd;cd;0.55;0.65;1};

[inputs,canceled]=inputsdlg(prompt,name,formats,defaultanswer);
if canceled==1; outstring=[]; return; end;
outstring={'[d2c c2d dmin dmax cmin cmax]=';['bleedthru_corr('...
    '''' inputs{1} '''' ',' '''' inputs{2} '''' ',' num2str(inputs{3}) ',' num2str(inputs{4}) ...
    num2str(inputs{5}) ');'];...
    'mask_out=mask_generator(output_files);';...
    'output_files=pearson_correction(d2c,c2d,cmin,cmax,dmin,dmax,mask_out,original_files);'};

end

