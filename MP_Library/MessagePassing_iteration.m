function [estFin, optFin] = ...
    MessagePassing_iteration(P,X,H_0,A,gOut, gU, gR, opt)

%% Setup
% Get options
nit     = opt.nit;              % number of iterations
nitMin  = opt.nitMin;           % minimum number of iterations
step    = opt.step;             % step size
stepMin = opt.stepMin;          % minimum step size
stepMax = opt.stepMax;          % maximum step size
adaptStep = opt.adaptStep;      % adaptive step size
stepIncr = opt.stepIncr;        % step inc on succesful step
stepDecr = opt.stepDecr;        % step dec on failed step
stepWindow = opt.stepWindow;    % step size check window size
tol = opt.tol;                  % Convergence tolerance
stepTol = opt.stepTol;          % minimum allowed step size
compVal = adaptStep;            % only compute cost function for adaptive
maxBadSteps = opt.maxBadSteps;  % maximum number of allowed bad steps
maxStepDecr = opt.maxStepDecr;  % amount to decrease maxStep after failures
WvarToPvarMax = opt.WvarToPvarMax;  % maximum Wvar/pvar ratio
%Get problem dimensions
A2=abs(A).^2;
X2=abs(X).^2;
P2=abs(P).^2;
%Assign Xvar and Svar mins
SvarMin = opt.SvarMin;

estFin.cov=false;
%% Initialization
state = [];
valIn = -inf;

rhat = opt.shat0;
rvar=opt.svar0;
uhat = opt.ghat0;
what = opt.what0;
uvar = opt.gvar0;
wvar=opt.wvar0;
zhat = opt.zhat0;
zvar=opt.zvar0;


%Placeholder initializations- values are not used
%Init valOpt empty
valOpt = [];
%Specify minimum variances
pvarMin = opt.pvarMin;
val = zeros(nit,1);
%% Iterations
% chat=what;
% cvar=wvar;
uhatBar=0;
whatBar=0;
rhatBar=0;
zhatBar=0;

%Control variable to end the iterations
stop = false;
it = 0;
failCount = 0;
step1 = 1;
shat = 0;
svar = 0;
zvarOpt = 0;
alphahat=0;
alphavar=0;
uhatBarOpt=0;
rhatBarOpt=0;
gammahat=0;
gammavar=0;
pvarOpt = 0;
zvar0Opt = 0;
uvarOpt=0;
muvarOpt=0;
betavarOpt=0;
% Main iteration loop
while ~stop
    % Iteration count
    it = it + 1;
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    betavar=zvar*X2;
    betavar=step1*betavar + (1-step1)*betavarOpt;
    betahat0=zhat*X;
    betahat=betahat0-gammahat.*betavar;
    if (compVal)
        valOut = sum(sum(gOut.logLike(betahat0,betavar)));
        val(it) = valOut + valIn;
    end
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
    if pass
        if any(isnan(betahat0(:))) || any(isinf(betahat0(:)))
            stop = true;
        else
            testVal2 = norm(rhat(:) - rhatBarOpt(:)) / norm(rhatBarOpt(:));
            testVal3 = norm(uhat(:) - uhatBarOpt(:)) / norm(uhatBarOpt(:));
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
        gammahatOpt = gammahat;
        gammavarOpt = gammavar;
        zhatBarOpt = zhatBar;
        zhatOpt = zhat;
        zvarOpt=zvar;
        betavarOpt=betavar;
        betavar = max(betavar, pvarMin);
        [qhat,qvar] = gOut.estim(betahat,betavar);
        pvarInv3 = 1 ./ betavar;
        gammahatNew = pvarInv3.*(qhat-betahat);
        gammavarNew = pvarInv3.*(1-min(qvar./betavar,WvarToPvarMax));
    end
    gammahat = (1-step1)*gammahatOpt + step1*gammahatNew;
    gammavar = (1-step1)*gammavarOpt + step1*gammavarNew;
    zhatBar = (1-step1)*zhatBarOpt + step1*zhatOpt;
    evar=1./(gammavar*X2');
    evar(evar > opt.varThresh) = opt.varThresh;
    ehat=zhatBar+evar.*(gammahat*X');
    evar = max(evar, SvarMin);
    
    
    what2 = abs(what).^2;
    uhat2 = abs(uhat).^2;
    zvar0= wvar*uhat2+what2*uvar;
    pvar=zvar0+wvar*uvar;
    
    pvar = step1*pvar + (1-step1)*pvarOpt;
    zvar0= step1*zvar0 + (1-step1)*zvar0Opt;
    zhat0 =what*uhat;
    phat = zhat0 - shat.*(zvar0);
    [zhat,zvar] =CAwgnEstimIn(phat,pvar).estim(ehat,evar);
    if pass
        pvarOpt = pvar;
        zvar0Opt = zvar0;
        pvar = max(pvar, pvarMin);
        uhatBarOpt = uhatBar;
        uhatOpt = uhat;
        whatBarOpt = whatBar;
        whatOpt = what;
        shatOpt = shat;
        svarOpt = svar;
        uvarOpt = uvar;
        wvarOpt = wvar;
        pvarInv = 1 ./ pvar;
        shatNew = pvarInv.*(zhat-phat);
        svarNew = pvarInv.*(1-min(zvar./pvar,WvarToPvarMax));
    end
    whatBar = (1-step1)*whatBarOpt + step1*whatOpt;
    shat = (1-step1)*shatOpt + step1*shatNew;
    svar = (1-step1)*svarOpt + step1*svarNew;
    uhatBar = (1-step1)*uhatBarOpt + step1*uhatOpt;
    bvar = 1./(abs(whatBar').^2*svar);
    bvar(bvar > opt.varThresh) = opt.varThresh;
    bGain = (1 - (bvar.*(wvar'*svar)));
    bGain = min(1,max(0,bGain));
    bhat = uhatBar.*bGain +bvar.*(whatBar'*shat);
    bvar = max(bvar, SvarMin);
    cvar = 1./(svar*(abs(uhatBar).^2)');
    cvar(cvar > opt.varThresh) = opt.varThresh;
    cGain = (1 - (cvar.*(svar*uvar')));
    cGain = min(1,max(0,cGain));
    chat = whatBar.*cGain +cvar.*(shat*uhatBar');
    cvar = max(cvar, SvarMin);
    [uhat,uvar,valInX] = gU.estim(bhat, bvar);
    muvar = A2*rvar*P2;
    muvar=step1*muvar + (1-step1)*muvarOpt;
    muhat=A*rhat*P-muvar.*alphahat;
    [what,wvar]=CAwgnEstimIn(muhat+H_0,muvar).estim(chat,cvar);
    if pass
        muvar = max(muvar, pvarMin);
        muvarOpt=muvar;
        rhatBarOpt = rhatBar;
        rhatOpt = rhat;
        rvarOpt = rvar;
        alphahatOpt = alphahat;
        alphavarOpt = alphavar;
        %         pvarInv2 = 1 ./ muvar;
        %         [nuhat,nuvar] = CAwgnEstimOut(chat-H_0, cvar).estim(muhat,muvar);
        %         alphahatNew = pvarInv2.*(nuhat-muhat);
        %         alphavarNew = pvarInv2.*(1-min(nuvar./muvar,WvarToPvarMax));
        alphavarNew=1./(cvar+muvar);
        alphahatNew=alphavarNew.*(chat-muhat-H_0);
    end
    alphahat = (1-step1)*alphahatOpt + step1*alphahatNew;
    alphavar = (1-step1)*alphavarOpt + step1*alphavarNew;
    rhatBar = (1-step1)*rhatBarOpt + step1*rhatOpt;
    dvar = 1./((A2)'*alphavar*(P2)');
    dvar(dvar > opt.varThresh) = opt.varThresh;
    dhat = rhatBar +dvar.*(A'*alphahat*P');
    dvar = max(dvar, SvarMin);
    [rhat,rvar,valInY] = gR.estim(dhat, dvar);
    uvar=step1*uvar + (1-step1)*uvarOpt;
    wvar=step1*wvar + (1-step1)*wvarOpt;
    rvar=step1*rvar + (1-step1)*rvarOpt;
    zvar=step1*zvar + (1-step1)*zvarOpt;
    valIn = sum ( valInX(:) )+sum( valInY(:) ) ;
    %Don't stop before minimum iteration count
    if it < nitMin
        stop = false;
    end
end
%% Save the final values
%Save the options object that was used
optFin = opt;
%Estimates of the two matrix factors
estFin.ghat = uhat;
estFin.gvar = uvar;
estFin.shat = rhat;
estFin.svar = rvar;
estFin.what=what;
estFin.wvar = wvar;
estFin.zhat = zhat;
estFin.zvar = zvar;
estFin.qhat = qhat;
estFin.qvar = qvar;

estFin.val=val(it);
