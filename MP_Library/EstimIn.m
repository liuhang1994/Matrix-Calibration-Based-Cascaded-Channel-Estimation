classdef EstimIn < hgsetget
    % EstimIn:  Base class for input function in the G-AMP algorithm
    %
    % Given input variables X, the function computes estimates
    % of X from observations of the form R = X + V, V ~ N(0,rvar)
    
    methods (Abstract)
        
        % Main estimation method:  For the sum-product algorithm,
        % the method should return:
        %
        %   xhat = E( X | R=rhat )
        %   xvar = var( X | R=rhat )
        %   val  = E( log ( p(X) / p(X|R=rhat) ) | R=rhat )
        %
        % For the max-sum algorithm:
        %
        %   xhat = argmax p(X|R=rhat)
        %   xvar = rvar./(1-fderiv2.*rvar)
        %   val =  log p(X) at X=xhat
        %
        % where fderiv2 = d^2/dx^2 log p(x) | x=xhat
        [xhat,xvar,val] = estim(obj,rhat,rvar)
       
        
        % Initial estimate:  Provides an initial estimate
        % Can be implemented in 
        [xhat,xvar,valInit] = estimInit(obj)
    end
    
    % Virtual methods that can be overloaded in the base class
    methods
        
        % Size operator.  Used for concatenation.  A value of nx=[]
        % implies that the size is auto-fit.  Autofit estimators cannot be
        % concatenated via the vertcat function
        function [nx,ncol] = size(obj)
            ncol = 1;
            nx = [];
        end
        
       % Vertical concatenation
        function obj =  vertcat(varargin)
            obj = EstimInConcat(varargin);
        end
        
    end
    
end