function [ outstring ] = Apply_tolerances_script()
%This is a test function for loading scripts into ForestFire. The code may
%look complex because I was testing a complex 'inputdlg' function, but all
%you need to do is just have your custom function output a string that will
%be executed by the cluster. 
%For example: outstring='disp(''hello world'');' would work just fine.
%   Philip Andresen, 8-22-2012 (august)

%function_apply_tolerances2_2color_batch(N_tol_min,N_tol_max,r0_tol_min,r0_tol_max,fr_un,max_unc,fname)

prompt={'Min. photons/molecule=','Max photons/molecule=',...
            'Min. PSF radius (r0)=','Max PSF radius (r0)=',...
            'Max uncertainty in PSF fit=',...
            'Max uncertainty in molecule position (nm)=',...
            'Input for this program is the previous output:'};
name='Apply Tolerances';
formats(1,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(2,1)=struct('type','edit','format','integer','limits',[1 inf],'style','edit');
formats(3,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(4,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(5,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(6,1)=struct('type','edit','format','float','limits',[0 inf],'style','edit');
formats(7,1)=struct('type','check','format','integer','limits',[0 1],'style','checkbox');
defaultanswer={150,4000,2.5,5,0.4,80,1};

[inputs,canceled]=inputsdlg(prompt,name,formats,defaultanswer);
if canceled==1; outstring=[]; return; end;
file='original_files';
if inputs{7}==1; file='output_files'; end;
outstring={['output_files = function_apply_tolerances2_2color_batch(' ...
    num2str(inputs{1}) ',' num2str(inputs{2}) ',' num2str(inputs{3}) ','...
    num2str(inputs{4}) ',' num2str(inputs{5}) ',' num2str(inputs{6}) ','...
    file ');']};

end

