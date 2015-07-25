function [ outstring ] = mask_generator_script()
%This is a test function for loading scripts into ForestFire. All
%you need to do is just have your custom function output a string that will
%be executed by the cluster. (A cell array of strings works too!)
%   Philip Andresen, 8-22-2012 (august)

%Special keywords you can use in any script?
%output_files
%original_files
%Generally the convention is to write your function like this:
%output_files=myfunction(input,input,input,original_files)
%Or write it like this if it depends on previous files:
%output_files=myfunction(input,input,input,output_files)
%
%The script system is really versatile, so you could even just write your
%code in ForestFire to look like this:
%A=myfunction(input,input,input,original_files)
%B=myfunction2(input,A);
%C=Anotherfunction(A,B);
%
%So long as your script function (of which this is an example) outputs a
%string or cell array of strings that can be executed as code, your
%functions will run just fine!


outstring={'mask_out=mask_generator(original_files)'};

end

