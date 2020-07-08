classdef EstimOut < hgsetget
% EstimOut:  Base class for output function in the GAMP algorithm
%
% For the sum-product algorithm, the output estimation function is defined 
% based on the output channel probability distribution p_{Y|Z}(y|z).  
% The value y is not passed to the function, but typically stored as a 
% member of the function.
%
% For the max-sum algorithm, the ouptut estimation function is based on an
% objective function fout(z).  For MAP estimation, this is generally set to
% the log likelihood:
% 
%   fout(z) = log P(y|z).    
    methods (Abstract)
        
        % Main estimation method:  For the sum-product algorithm,
        % the method should return:
        %
        %   zhat = E( Z | Y)
        %   zvar = var( Z | Y )
        %
        % where Z = N(phat, pvar).  For the max-sum algorithm, the method
        % should return:
        %   
        %   zhat = argmax_z [ fout(z) - (1/2*pvar)*abs( z-p )^2 ]
        %   zvar = pvar ./( 1 - fout''(zhat)*pvar) )
        %
        % Note that if fout is concave, fout''(zhat) <= 0, so 
        %   0 <= zvar <= pvar.
        [zhat,zvar] = estim(obj,phat,pvar)
        
        % Log-likelihood:  For sum-product, the method should return 
        %   E( log p_{Y|Z}(y|Z) )  with Z = N(zhat,zvar)
        %
        % For max-sum, it should return fout(zhat).  The variable zvar is
        % ignored.
        ll = logLike(obj,zhat,zvar) 
        
        % Compute output cost:
        % For sum-product (real case) compute
        %   (Axhat-phatfix)^2/(2*pvar) + log int_z p_{Y|Z}(y|z) N(z;phatfix, pvar) 
        %   with phatfix such that Axhat=estim(phatfix,pvar).
        % For sum-product (complex case) compute
        %   abs(Axhat-phatfix)^2/pvar + log int_z p_{Y|Z}(y|z) CN(z;phatfix, pvar) 
        % For max-sum GAMP, compute
        %   log p_{Y|Z}(y|z) @ z = Axhat
        scale = logScale(obj,Axhat,pvar,phat)        
        
    end
    
    % Virtual functions that may be overwritten if desired
    methods
        
        % Return number of columns
        function S = numColumns(obj) %#ok<MANU>
            
            %By default, return 1
            S = 1;
        end
        
        % Size.  This needs to be overwritten to support vertcat.
        function [nz,ncol] = size(obj)
            nz = [];
            ncol = 1;            
        end
        
        % Vertical concatenation
        function obj =  vertcat(varargin)
            obj = EstimOutConcat(varargin);
        end
        
    end
    
end