classdef MPOpt 
% structure that contains the algorithm control parameters
    properties
        
        
        %Maximum number of iterations
        nit = 1500;
        %Minimum number of iterations
        nitMin = 0; %0 for no effect
        
        %Convergence tolerance.
        tol = 1e-4;
        
 
        
        
        
        %***** Initialization
        %Provide initializations for the inital variables
        shat0 = [];
        ghat0 = [];
        svar0 = [];
        gvar0 = [];
        what0 =[];
        wvar0=[];
        zhat0=[];
        zvar0=[];
        
        
        %***** Step Size Control
        %Logical flag to include a step size in the pvar/Wvar calculation.
        %This momentum term often improves numerical performance. On by
        %defualt.
        pvarStep = true;
        %Initial step size, or fixed size for non-adaptive steps
        step = 0.05;
        % Enable adaptive step size
        adaptStep = true;
        % Minimum step size
        stepMin = 0.05;
        % Maximum step size
        stepMax = 0.5;
        % Multiplicative step size increase, when successful
        stepIncr = 1.1;
        % Multiplicative step size decrease, when unsuccessful
        stepDecr = 0.5;
        %Maximum number of allowed failed steps before we decrease stepMax,
        %inf for no effect
        maxBadSteps = inf;
        %Amount to decrease stepMax after maxBadSteps failed steps, 1 for
        %no effect
        maxStepDecr = 0.5;
        
        %Create a window for the adaptive step size test. Setting this to
        %zero causes it to have no effect. For postive integer values,
        %creats a moving window of this length when checking the step size
        %acceptance criteria. The new value is only required to be better
        %than the worst in this window, i.e. the step size criteria is not
        %required to monotonicaly increase.
        stepWindow = 10;
        
        % Iterations are termined when the step size becomes smaller
        % than this value. Set to -1 to disable
        stepTol = -1;
        
        %This is a filter value on the steps. It slows down the early steps
        %in a smooth way. Set to less than 1 for no effect, In particular,
        %the steps are computed as step1 = step*it/(it + stepFilter)
        stepFilter = 0;
        
        
        %***** Variance Control
        %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        
        %Variance lower bound. reset to this lower bound if the variance terms
        %are too less
        pvarMin = 1e-5;
        SvarMin = 1e-5;
        WvarToPvarMax = 0.99;   % prevents negative variance, should be near 1
        
        %Variance upper bound. Prevent variances to be too large
        varThresh = 1e6; 
    end
    
    methods
        
        % Constructor with default options
        function opt = MPOpt(varargin)
            if nargin == 0
                % No custom parameters values, thus create default object
                return
            elseif mod(nargin, 2) == 0
                % User is providing property/value pairs
                names = fieldnames(opt);    % Get names of class properties
                
                % Iterate through every property/value pair, assigning
                % user-specified values.  Note that the matching is NOT
                % case-sensitive
                for i = 1 : 2 : nargin - 1
                    if any(strcmpi(varargin{i}, names))
                        % User has specified a valid property
                        propName = names(strcmpi(varargin{i}, names));
                        opt.(propName{1}) = varargin{i+1};
                    else
                        % Unrecognized property name
                        error('BDMPOpt: %s is an unrecognized option', ...
                            num2str(varargin{i}));
                    end
                end
                return
            else
                error(['The BDMPOpt constructor requires arguments ' ...
                    'be provided in pairs, e.g., BDMPOpt(''adaptStep'',' ...
                    ' false, ''nit'', 50)'])
            end
        end
    end
    
end
