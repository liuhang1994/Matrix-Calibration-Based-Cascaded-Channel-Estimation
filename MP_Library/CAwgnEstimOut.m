classdef CAwgnEstimOut < EstimOut
    % This script is from GAMPMATLAB package.
    
    % CAwgnEstimOut:  CAWGN scalar output estimation function
    %
    % Corresponds to an output channel of the form
    %   y = scale*z + CN(0, wvar)
    
    properties
        y;                 % Measured output
        wvar;              % Variance
        scale = 1;         % scale factor
        maxSumVal = false; % True indicates to compute output for max-sum
        autoTune = false;   % Set to true for tuning of mean and/or variance
        disableTune = false;% Set to true to temporarily disable tuning
        tuneMethod = 'Bethe';  % Tuning method, in {ML,Bethe,EM0,EM}
        tuneDim = 'joint';  % Dimension to autoTune over, in {joint,col,row}
        tuneDamp = 0.1;     % Damping factor for autoTune in (0,1]
        counter = 0;        % Counter to delay tuning
        wvar_min = 1e-20;   % Minimum allowed value of wvar    
    end
    
    methods
        % Constructor
        function obj = CAwgnEstimOut(y, wvar, maxSumVal, varargin)
            obj = obj@EstimOut;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.y = y;
                obj.wvar = wvar;
                if (nargin >= 3)
                    if (~isempty(maxSumVal))
                        obj.maxSumVal = maxSumVal;
                    end
                end
                if (nargin >= 4)
                    if isnumeric(varargin{1})
                        % make backwards compatible: 4th argument can specify the scale
                        obj.scale = varargin{1};
                    else
                        for i = 1:2:length(varargin)
                            obj.(varargin{i}) = varargin{i+1};
                        end
                    end
                end
                
                %Warn user about zero-valued noise variance
                %if any(obj.wvar == 0)
                %    warning(['Tiny non-zero variances will be used for'...
                %        ' computing log likelihoods. May cause problems'...
                %        ' with adaptive step size if used.']) %#ok<*WNTAG>
                %end
            end
        end

        % Set methods
        function obj = set.y(obj, y)
            obj.y = y;
        end

        function obj = set.wvar(obj, wvar)
            assert(all(wvar(:) >= 0), ...
                'CAwgnEstimOut: wvar must be non-negative');
            obj.wvar = wvar;
        end

        function obj = set.wvar_min(obj, wvar_min)
            assert(all(wvar_min(:) > 0), ...
                'CAwgnEstimOut: wvar_min must be positive');
            obj.wvar_min = wvar_min;
        end

        function obj = set.maxSumVal(obj, maxSumVal)
            assert(isscalar(maxSumVal)&&(ismember(maxSumVal,[0,1])||islogical(maxSumVal)), ...
                'CAwgnEstimOut: maxSumVal must be a logical scalar');
            obj.maxSumVal = maxSumVal;
        end

        function obj = set.scale(obj, scale)
            assert(isnumeric(scale)&&isscalar(scale)&&(scale>0), ...
                'CAwgnEstimOut: scale must be a positive scalar');
            obj.scale = scale;
        end

        function set.disableTune(obj, flag)
            assert(isscalar(flag)&&(ismember(flag,[0,1])||islogical(flag)), ...
                'CAwgnEstimOut: disableTune must be a logical scalar');
            obj.disableTune = flag;
        end
        
        % Size
        function [nz,ncol] = size(obj)
            [nz,ncol] = size(obj.y);
        end
        
        % AWGN estimation function
        % Provides the posterior mean and variance of variable z
        % from an observation y = scale*z + w 
        % where z = CN(phat,pvar), w = CN(0,wvar)
        function [zhat, zvar, partition] = estim(obj, phat, pvar)
            
            % Extract quantities
            y = obj.y;
            scale = obj.scale;
            scale2pvar = (scale^2)*pvar;

            % Compute posterior mean and variance
            wvar = obj.wvar;
            gain = pvar./(scale2pvar + wvar);
            zhat = (scale*gain).*(y-scale*phat) + phat;
            zvar = wvar.*gain;

            % Compute partition function
            if nargout==3
                partition = (1./(pi*(scale2pvar+wvar))).*exp( ...
                                -(abs(phat-y).^2) ./ (scale2pvar+wvar) );
            end
            
            % Tune noise variance
            if obj.autoTune && ~obj.disableTune
                if (obj.counter>0), % don't tune yet
                    obj.counter = obj.counter-1; % decrement counter 
                else % tune now
                    [M,T] = size(phat);
                    damp = obj.tuneDamp;
                    %Learn variance, averaging over columns and/or rows
                    switch obj.tuneMethod
                      case 'ML' % argmax_wvar Z(y;wvar)=\int p(y|z;wvar) N(z;phat,pvar) dz
                        switch obj.tuneDim
                          case 'joint'
                            wvar1 = mean(abs(scale*phat(:)-y(:)).^2 - scale2pvar(:));
                            wvar0 = max(obj.wvar_min, wvar1)*ones(size(obj.wvar));
                          case 'col'
                            wvar1 = (1/M)*sum(abs(scale*phat-y).^2 - scale2pvar,1);
                            wvar0 = ones(M,1)*max(obj.wvar_min, wvar1);
                          case 'row'
                            wvar1 = (1/T)*sum(abs(scale*phat-y).^2 - scale2pvar,2);
                            wvar0 = max(obj.wvar_min, wvar1)*ones(1,T);
                          otherwise
                            error('Invalid tuning dimension in CAwgnEstimOut');
                        end
                        if damp==1
                          obj.wvar = wvar0;
                        else % apply damping
                          obj.wvar = exp( (1-damp)*log(obj.wvar) + damp*log(wvar0));
                        end

                      case 'Bethe' % Method from Krzakala et al J.Stat.Mech. 2012
                        svar = 1./(scale2pvar + obj.wvar);
                        shat = (y-scale*phat).*svar;
                        switch obj.tuneDim
                          case 'joint'
                            ratio = sum(abs(shat(:)).^2)/sum(svar(:));
                            if damp~=1, ratio = ratio.^damp; end;
                            obj.wvar = max(obj.wvar_min, obj.wvar*ratio);
                          case 'col'
                            ratio = sum(abs(shat).^2,1)./sum(svar,1);
                            if damp~=1, ratio = ratio.^damp; end;
                            obj.wvar = max(obj.wvar_min, obj.wvar.*(ones(M,1)*ratio));
                          case 'row'
                            ratio = sum(abs(shat).^2,2)./sum(svar,2);
                            if damp~=1, ratio = ratio.^damp; end;
                            obj.wvar = max(obj.wvar_min, obj.wvar.*(ratio*ones(1,T)));
                          otherwise
                            error('Invalid tuning dimension in CAwgnEstimOut');
                        end

                      case 'EM0'
                        switch obj.tuneDim
                          case 'joint'
                            wvar1 = mean(abs(y(:)-zhat(:)).^2);
                            wvar0 = max(obj.wvar_min, wvar1)*ones(size(obj.wvar));
                          case 'col'
                            wvar1 = (1/M)*sum(abs(y-zhat).^2,1);
                            wvar0 = ones(M,1)*max(obj.wvar_min, wvar1);
                          case 'row'
                            wvar1 = (1/T)*sum(abs(y-zhat).^2,2);
                            wvar0 = max(obj.wvar_min, wvar1)*ones(1,T);
                          otherwise
                            error('Invalid tuning dimension in CAwgnEstimOut');
                        end
                        if damp==1
                          obj.wvar = wvar0;
                        else % apply damping
                          obj.wvar = exp( (1-damp)*log(obj.wvar) + damp*log(wvar0));
                        end

                      case 'EM'
                        switch obj.tuneDim
                          case 'joint'
                            wvar1 = mean(abs(y(:)-zhat(:)).^2 + zvar(:));
                            wvar0 = max(obj.wvar_min, wvar1)*ones(size(obj.wvar));
                          case 'col'
                            wvar1 = (1/M)*sum(abs(y-zhat).^2 + zvar,1);
                            wvar0 = ones(M,1)*max(obj.wvar_min, wvar1);
                          case 'row'
                            wvar1 = (1/T)*sum(abs(y-zhat).^2 + zvar,2);
                            wvar0 = max(obj.wvar_min, wvar1)*ones(1,T);
                          otherwise
                            error('Invalid tuning dimension in CAwgnEstimOut');
                        end
                        if damp==1
                          obj.wvar = wvar0;
                        else % apply damping
                          obj.wvar = exp( (1-damp)*log(obj.wvar) + damp*log(wvar0));
                        end

                      otherwise
                        error('Invalid tuning method in CAwgnEstimOut');
                    end
                end
            end

        end
        
        % Compute log likelihood
        % For sum-product GAMP, compute
        %   E( log p_{Y|Z}(y|z) ) with z = CN(phat, pvar)
        % For max-sum GAMP compute
        %   log p_{Y|Z}(y|z) @ z = phat
        function ll = logLike(obj,phat,pvar)
            
            % Ensure variance is small positive number
            wvar_pos = max(obj.wvar_min, obj.wvar);
            
            % Get scale
            scale = obj.scale;            

            % Compute log-likelihood
            if ~(obj.maxSumVal)

                predErr = (abs(obj.y-phat).^2 + pvar)./wvar_pos;
            else
                predErr = (abs(obj.y-scale*phat).^2)./wvar_pos;
            end
            ll = -(predErr); %return the values without summing
        end
        
        % Compute output cost:
        % For sum-product compute
        %   abs(Axhat-phatfix)^2/(pvar) + log int_z p_{Y|Z}(y|z) CN(z;phatfix, pvar) 
        %   with phatfix such that Axhat=estim(phatfix,pvar).
        % For max-sum GAMP, compute
        %   log p_{Y|Z}(y|z) @ z = Axhat
        function ll = logScale(obj,Axhat,pvar,phat)
                   
            % Ensure variance is small positive number
            wvar1 = max(obj.wvar_min, obj.wvar);
            
            %Get the scale
            s = obj.scale;   
           
            % Compute output cost
            if ~(obj.maxSumVal)

                % Compute output cost
                closed_form = true;
                if closed_form
                    
                    %Closed form update
                    ll = -log(abs(s)^2*pvar + wvar1) ...
                        - abs(obj.y - s*Axhat).^2./wvar1 - log(pi); 
                else
                    % Find the fixed-point of phat
                    opt.phat0 = Axhat; % works better than phat
                    opt.alg = 1; % approximate newton's method
                    opt.maxIter = 3; 
                    opt.tol = 1e-4; 
                    opt.stepsize = 1; 
                    opt.regularization = obj.wvar^2;  % works well up to SNR=160dB
                    opt.debug = false;
                    phatfix = estimInvert(obj,Axhat,pvar,opt);

                    % Compute log int_z p_{Y|Z}(y|z) CN(z;phatfix, pvar)
                    ls = -log(pi*(obj.wvar + abs(s)^2*pvar)) ...
                        - abs(obj.y - s*phatfix).^2 ./ (obj.wvar + abs(s)^2*pvar);

                    % Combine to form output cost
                    ll = ls + abs(Axhat - phatfix).^2./pvar;
                end;

            else
                % Output cost is simply the log likelihood
                ll = -abs(obj.y-s*Axhat).^2./wvar1;
            end 
            
        end
        
        function S = numColumns(obj)
            %Return number of columns of Y
            S = size(obj.y,2);
        end
        
        % Generate random samples from p(y|z)
        function y = genRand(obj, z)
            y = sqrt(obj.wvar/2).*randn(size(z)) + ...
                1j*sqrt(obj.wvar/2).*randn(size(z)) + obj.scale.*z;
        end
    end
    
end

