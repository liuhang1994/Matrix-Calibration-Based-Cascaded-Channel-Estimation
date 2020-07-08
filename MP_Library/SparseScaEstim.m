classdef SparseScaEstim < EstimIn
    % SparseScaEstim:  Scalar estimator class with sparsity
    %
    % Let baseEstim be the estimator for a base random variable X1.
    % Then SparseScaEstim is the estimator for a new random variable:
    %   X = X1 with prob p1
    %     = x0 with prob 1-p1
    % where x0 is a constant (default = 0)
    properties
        x0 = 0;             % location of Dirac delta
        p1 = 0.5;           % Sparsity rate
        estim1;             % Base estimator when U=1
        autoTune = false;   % tuning of sparsity parameter?
        disableTune = false;% temporarily disable tuning of this and components?
        tuneDim = 'joint';  % Determine dimension to autoTune over, must be 
                            % 'joint', 'col', or 'row'
        counter = 0;        % Counter to delay tuning
    end
    
    properties (Hidden)
        LogLikeFlag = false; % Indicates if estim1 has a loglikey method
        weightFlag = false;  % Indicates if estim1 has a weightFlag property
    end
    
    methods
        % Constructor
        function obj = SparseScaEstim(estim1, p1, x0, varargin)
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.p1 = p1;
                obj.estim1 = estim1;
                if (nargin >= 3 && ~isempty(x0))
                    obj.x0 = x0;
                end
            end
            for i = 1:2:length(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end
        
        % Set method for estim1
        function set.estim1(obj, Estim1)
            % Check to ensure input Estim1 is a valid EstimIn class
            if isa(Estim1, 'EstimIn')
                obj.LogLikeFlag = ismethod(Estim1, 'loglikey'); %#ok<MCSUP>
                obj.weightFlag = ~isempty(findprop(Estim1, 'mixWeight'));
                obj.estim1 = Estim1;
            else
                error('estim1 must be a valid EstimIn object')
            end
        end
        
        % Set method for disableTune
        function set.disableTune(obj, flag)
            if islogical(flag)
                % change property of this class 
                obj.disableTune = flag;
                % change property of component class 
                if any(strcmp('disableTune', properties(obj.estim1)))
                    obj.estim1.disableTune = flag;
                end
            else
                error('disableTune must be a boolean')
            end
        end
        
        % Compute prior mean and variance
        function [xhat, xvar, valInit] = estimInit(obj)
            [xhat1, xvar1, valInit1] = obj.estim1.estimInit;
            xhat = obj.p1.*xhat1 + (1-obj.p1).*obj.x0;
            xvar = obj.p1.*(abs(xhat1).^2 + xvar1) ...
                + (1-obj.p1).*(abs(obj.x0).^2) - abs(xhat).^2;
            valInit = obj.p1.*valInit1;
        end
        
        % Compute posterior mean and variance from Gaussian estimate
        function [xhat, xvar, klDivNeg, py1] = estim(obj, rhat, rvar)
            
            % Compute the activity probabilities
            if ~obj.LogLikeFlag
                % This EstimIn class does not implement the loglikey
                % method, thus convert from probability domain to log-prob
                % domain
                
                % Get log-likelihood of rhat for U=1 and U=0
                loglike1 = log( obj.estim1.plikey(rhat, rvar) );

            else
                % This EstimIn class implements the loglikey method, thus
                % work in the log-domain
                
                % Get log-likelihood of rhat for U=1 and U=0
                loglike1 = obj.estim1.loglikey(rhat, rvar);
            end
            
            %Handle real and complex cases separately
            rvar(rvar < eps) = eps;     % for numerical stability...
            if isreal(rhat)
                loglike0 = -0.5*( log(2*pi) + log(rvar) + ...
                    ((rhat - obj.x0).^2)./rvar );
            else
                loglike0 = -( log(pi) + log(rvar) + ...
                    (abs(rhat - obj.x0).^2)./rvar );
            end
            
            % Convert log-domain quantities into posterior activity
            % probabilities (i.e., py1 = Pr{X=X1 | y}, py0 = Pr{X=x0 | y})
            exparg = loglike0 - loglike1 + log(1 - obj.p1) - log(obj.p1);
            maxarg = 500; 
            exparg = max(min(exparg,maxarg),-maxarg); % numerical robustness
            py1 = (1 + exp(exparg)).^(-1);
            py0 = 1 - py1;
            
            %Update the sparsity rate parameter
            if obj.autoTune && ~obj.disableTune

              if (obj.counter>0), % don't tune yet
                obj.counter = obj.counter-1; % decrement counter 
              else % tune now

                [N, T] = size(rhat);
                %Average over all elements, per column, or per row
                p1_old = obj.p1;
                switch obj.tuneDim
                    case 'joint'
                        obj.p1 = sum(py1(:))/N/T;
                    case 'col'
                        obj.p1 = repmat(sum(py1)/N,[N 1]);
                    case 'row'
                        obj.p1 = repmat(sum(py1,2)/T, [1 T]);
                    otherwise 
                        error('Invalid tuning dimension in SparseScaEstim');
                end
                if obj.weightFlag
                    set(obj.estim1, 'mixWeight', py1./obj.p1)

                    %Try keeping p1*var0 constant; assumes that estim1 has var0!
                    if 0
                      set(obj.estim1, 'var0', obj.estim1.var0.*p1_old./obj.p1) 
                      %p1var0 = obj.p1*obj.estim1.var0 
                    end

                else
                    if ~isempty(findprop(obj.estim1,'autoTune'))
                      if obj.estim1.autoTune
                        warning(strcat('Auto tuning may fail: ',...
                                class(obj.estim1),...
                                ' does not include property mixWeight.'))
                      end 
                    end
                end

              end
            end  
            
            % Compute mean and variance
            if (nargout >= 3)
                [xhat1, xvar1, klDivNeg1] = obj.estim1.estim(rhat, rvar);
            else
                [xhat1, xvar1] = obj.estim1.estim(rhat, rvar);
            end
            xhat = py1.*xhat1 + py0.*obj.x0;
%           xvar = py1.*(abs(xhat1).^2 + xvar1) + py0.*(abs(obj.x0).^2)...
%               - abs(xhat).^2;
            xvar = py1.*(abs(xhat1).^2 - abs(xhat).^2) + py1.*xvar1 ... 
                + py0.*(abs(obj.x0).^2 -abs(xhat).^2);
            
            % Compute negative K-L divergence
            if (nargout >= 3)
                klDivNeg = py1.*klDivNeg1 ...
                    + py1.*log(max(1e-8,obj.p1)./max(py1,1e-8)) ...
                    + py0.*log(max(1e-8,(1-obj.p1))./max(py0,1e-8));
            end
            
        end
        
        % Generate random samples
        function x = genRand(obj, nx)
            x1 = obj.estim1.genRand(nx);
            p = rand(size(x1)) < obj.p1;
            x = x1.*p + (1-p).*obj.x0;
        end
        
        % Get the points in the distribution
        function x0 = getPoints(obj)
            x0 = [0; obj.estim1.getPoints()];
        end
        
        % Set sparsity level
        function setSparseProb(obj, p1)
            obj.p1 = p1;
        end
        
    end
    
end

