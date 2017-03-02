classdef parameter < handle
    %CLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( SetAccess = public)
        T %table of parameters
    end
    
    methods
      function obj = parameter(table)
          obj.T = table;
      end
    end
    
end

