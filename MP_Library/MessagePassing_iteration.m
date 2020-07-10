function [estFin, optFin] = ...
    MessagePassing_iteration(P,X,H_0,A,gOut, gU, gR, opt)
% Message passing iterations in Algorithm 1
% Quantities are computed in a matrix form
% Setup

% Get options
nit     = opt.nit;              % number of iterations
nitMin  = opt.nitMin;           % minimum number of iterations
tol = opt.tol;                  % Convergence tolerance
WvarToPvarMax = opt.WvarToPvarMax;  % maximum Wvar/pvar ratio

% use for adaptive damping, details please refer to Ref. [23]
step    = opt.step;             % step size
stepMin = opt.stepMin;          % minimum step size
stepMax = opt.stepMax;          % maximum step size
adaptStep = opt.adaptStep;      % adaptive step size
stepIncr = opt.stepIncr;        % step inc on succesful step
stepDecr = opt.stepDecr;        % step dec on failed step
stepWindow = opt.stepWindow;    % step size check window size
stepTol = opt.stepTol;          % minimum allowed step size
compVal = adaptStep;            % only compute cost function for adaptive
maxBadSteps = opt.maxBadSteps;  % maximum number of allowed bad steps
maxStepDecr = opt.maxStepDecr;  % amount to decrease maxStep after failures

%Compute the magnitude square matrices \abs{*}.^2
A2=abs(A).^2;
X2=abs(X).^2;
P2=abs(P).^2;
%Specify minimum variances
pvarMin = opt.pvarMin;
SvarMin = opt.SvarMin;
%whether early stop occurs
estFin.cov=false;
% Initialization
state = [];
valIn = -inf;
shat=  opt.shat0;
svar=  opt.svar0;
ghat=  opt.ghat0;
what=  opt.what0;
gvar=  opt.gvar0;
wvar=  opt.wvar0;
zhat=  opt.zhat0;
zvar=  opt.zvar0;
valOpt = [];
val = zeros(nit,1);




% Inital intermediate variables
ghatBar=0;
whatBar=0;
zhatBar=0;
zvarOpt = 0;
gvarOpt=0;
ghatBarOpt=0;



gammahat=0;
gammavar=0;
xihat = 0;
xivar = 0;


pvarOpt = 0;
zvar0Opt = 0;
betavarOpt=0;
muvarOpt=0;



shatBar=0;
shatBarOpt=0;
alphahat=0;
alphavar=0;

%Control variable to end the iterations
stop = false;
it = 0;
failCount = 0;

% inital step size
step1 = 1;
% Main iteration loop
while ~stop
    
    % Iteration count
    it = it + 1;
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    
    betavar=zvar*X2;  %v^beta, eq.(50a)
    betavar=step1*betavar + (1-step1)*betavarOpt;  %damp betavar
    betahat0=zhat*X;  %\hat beta without onsager correction
    betahat=betahat0-gammahat.*betavar; %\hat beta, eq. (50b)
    
    
    % adaptive damping: compute the output objective
    if (compVal)
        valOut = sum(sum(gOut.logLike(betahat0,betavar)));
        val(it) = valOut + valIn;
    end
    
    % check whether this iteration is good enough to pass
    if ~isempty(valOpt)
        stopInd = length(valOpt);
        startInd = max(1,stopInd - stepWindow);
        pass = (val(it) > min(valOpt(startInd:stopInd))) ||...
            ~adaptStep || (step <= stepMin);
    else
        pass = true;
    end
    
    if pass
        %Slightly inrease step size after pass if using adaptive steps
        step = stepIncr*step;
        valOpt = [valOpt val(it)];
        step = min([max([step stepMin]) stepMax]);
    else
        failCount = failCount + 1;
        if failCount > maxBadSteps
            failCount = 0;
            stepMax = max(stepMin,maxStepDecr*stepMax);
        end
        step = max(stepMin, stepDecr*step);
        if step < stepTol
            stop = true;
        end
    end
    % check stop condition
    if pass
        if any(isnan(betahat0(:))) || any(isinf(betahat0(:)))
            stop = true;
        else
            testVal2 = norm(shat(:) - shatBarOpt(:)) / norm(shatBarOpt(:));
            testVal3 = norm(ghat(:) - ghatBarOpt(:)) / norm(ghatBarOpt(:));
            if (it > 1) &&  (testVal3< tol) &&(testVal2<tol)
                stop = true;
                estFin.cov=true;
            end
        end
    end
    
    if it > 1 || ~isempty(state)
        step1 = step;
    end
    
    
    if pass
        
        % *Opt means the quantity in the previous iteration, used for
        % damping
        gammahatOpt = gammahat;
        gammavarOpt = gammavar;
        zhatBarOpt = zhatBar;
        zhatOpt = zhat;
        zvarOpt=zvar;
        betavarOpt=betavar;
        betavar = max(betavar, pvarMin);
        
        % qvar= tau_N*betavar/(betavar+tau_N)
        % qhat= (y*betavar+tau_N*betahat)/(betavar+tau_N)
        % two intermediate variables
        [qhat,qvar] = gOut.estim(betahat,betavar);
        pvarInv3 = 1 ./ betavar;
        gammahatNew = pvarInv3.*(qhat-betahat); % \hat gamma, eq.(49a)
        gammavarNew = pvarInv3.*(1-min(qvar./betavar,WvarToPvarMax)); % v^\gamma, eq. (49b)
    end
    %damp gammahat and gammavar
    gammahat = (1-step1)*gammahatOpt + step1*gammahatNew;
    gammavar = (1-step1)*gammavarOpt + step1*gammavarNew;
    %damp zhat
    zhatBar = (1-step1)*zhatBarOpt + step1*zhatOpt;
    
    
    evar=1./(gammavar*X2'); %v^e, eq. (49c)
    evar(evar > opt.varThresh) = opt.varThresh;
    ehat=zhatBar+evar.*(gammahat*X'); %\hat e, eq. (49d)
    evar = max(evar, SvarMin);
    
    
    %magnitude square 
    what2 = abs(what).^2;
    ghat2 = abs(ghat).^2;
    
    % v^p, eq. (56a)
    zvar0= wvar*ghat2+what2*gvar;
    pvar=zvar0+wvar*gvar;
    % damp pvar
    zvar0= step1*zvar0 + (1-step1)*zvar0Opt;
    pvar = step1*pvar + (1-step1)*pvarOpt;
    
    % \hat p, eq. (56b)
    zhat0 =what*ghat;
    phat = zhat0 - xihat.*(zvar0);
    
    % v^z and \hat z, eqs. (48c)--(48d)
    [zhat,zvar] =CAwgnEstimIn(phat,pvar).estim(ehat,evar);
    
    if pass
        % *Opt means the quantity in the previous iteration, used for
        % damping
        pvarOpt = pvar;
        zvar0Opt = zvar0;
        pvar = max(pvar, pvarMin);
        whatBarOpt = whatBar;
        whatOpt = what;
        wvarOpt = wvar;
        
        
        
        ghatBarOpt = ghatBar;
        ghatOpt = ghat;
        gvarOpt = gvar;      
        xihatOpt = xihat;
        xivarOpt = xivar;
        %\hat xi and v^\xi, eq. (54)
        pvarInv = 1 ./ pvar;
        xihatNew = pvarInv.*(zhat-phat);
        xivarNew = pvarInv.*(1-min(zvar./pvar,WvarToPvarMax));
    end
    %damp xihat and xivar
    xihat = (1-step1)*xihatOpt + step1*xihatNew;
    xivar = (1-step1)*xivarOpt + step1*xivarNew;
    %damp what and wvar
    whatBar = (1-step1)*whatBarOpt + step1*whatOpt;
    %damp ghat and gvar
    ghatBar = (1-step1)*ghatBarOpt + step1*ghatOpt;
    
    %eqs. (53a)--(53b)
    bvar = 1./(abs(whatBar').^2*xivar);
    bvar(bvar > opt.varThresh) = opt.varThresh;
    bGain = (1 - (bvar.*(wvar'*xivar)));
    bGain = min(1,max(0,bGain));
    bhat = ghatBar.*bGain +bvar.*(whatBar'*xihat);
    bvar = max(bvar, SvarMin);
    
    %eqs. (53c)--(53d)
    cvar = 1./(xivar*(abs(ghatBar).^2)');
    cvar(cvar > opt.varThresh) = opt.varThresh;
    cGain = (1 - (cvar.*(xivar*gvar')));
    cGain = min(1,max(0,cGain));
    chat = whatBar.*cGain +cvar.*(xihat*ghatBar');
    cvar = max(cvar, SvarMin);
    
    
    %eq. (55)
    [ghat,gvar,valInX] = gU.estim(bhat, bvar);
    
    
    %eqs. (61a)--(61b)
    muvar = A2*svar*P2;
    muvar=step1*muvar + (1-step1)*muvarOpt;
    muhat=A*shat*P-muvar.*alphahat;
    %eq. (63)
    [what,wvar]=CAwgnEstimIn(muhat+H_0,muvar).estim(chat,cvar);
    if pass
        muvarOpt=muvar;
        muvar = max(muvar, pvarMin);
        
        shatBarOpt = shatBar;
        shatOpt = shat;
        svarOpt = svar;
        alphahatOpt = alphahat;
        alphavarOpt = alphavar;
        
        
        
        %         pvarInv2 = 1 ./ muvar;
        %         [nuhat,nuvar] = CAwgnEstimOut(chat-H_0, cvar).estim(muhat,muvar);
        %         alphahatNew = pvarInv2.*(nuhat-muhat);
        %         alphavarNew = pvarInv2.*(1-min(nuvar./muvar,WvarToPvarMax));
        
        
        %eq. (62)
        alphavarNew=1./(cvar+muvar);
        alphahatNew=alphavarNew.*(chat-muhat-H_0);
    end
    %damping
    alphahat = (1-step1)*alphahatOpt + step1*alphahatNew;
    alphavar = (1-step1)*alphavarOpt + step1*alphavarNew;
    shatBar = (1-step1)*shatBarOpt + step1*shatOpt;
    
    %eqs. (61c)--(61d)
    dvar = 1./((A2)'*alphavar*(P2)');
    dvar(dvar > opt.varThresh) = opt.varThresh;
    dhat = shatBar +dvar.*(A'*alphahat*P');
    dvar = max(dvar, SvarMin);
    
    %eq. (64)
    [shat,svar,valInY] = gR.estim(dhat, dvar);
    
    %damp variance terms
    gvar=step1*gvar + (1-step1)*gvarOpt;
    wvar=step1*wvar + (1-step1)*wvarOpt;
    svar=step1*svar + (1-step1)*svarOpt;
    zvar=step1*zvar + (1-step1)*zvarOpt;
    valIn = sum ( valInX(:) )+sum( valInY(:) ) ;
    %Don't stop before minimum iteration count
    if it < nitMin
        stop = false;
    end
end


% Save the final values
%Save the options object that was used
optFin = opt;
%Estimates of the two matrix factors
estFin.ghat = ghat;
estFin.gvar = gvar;
estFin.shat = shat;
estFin.svar = svar;
estFin.what=what;
estFin.wvar = wvar;
estFin.zhat = zhat;
estFin.zvar = zvar;
estFin.qhat = qhat;
estFin.qvar = qvar;

estFin.val=val(it);
