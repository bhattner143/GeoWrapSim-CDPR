% Class to store the configuration of different robots from the XML files
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description    :
%    This class 

classdef CWOptWithGenIntDetBase < handle
    %UNTITLED17 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        wrap_model_config
    end
    methods
        %%
        function obj = CWOptWithGenIntDetBase(wrap_model_config)
            %UNTITLED17 Construct an instance of this class
            obj. wrap_model_config =  wrap_model_config;
        end
    end
end