function [MSE_S,MSE_G]=replica_iteration(snr,...
    K,M,M_prime,T,L,L_prime,lambda_S,lambda_G,tau_S,tau_G,tau_X,tau_H0)
% iterative algorithm for replica analysis
int_region=15; % all the numerical intergrations are over [-int_region,int_region]
threshold=1e-5; % early stop if the change is less than threshold
maxite=1000; % maximum number of iterations

MSE_S=zeros(1,length(snr));
MSE_G=zeros(1,length(snr));
for i=1:length(snr)
    % eq. (34)
    Q_S=lambda_S*tau_S;
    Q_G=lambda_G*tau_G;
    Q_W=L_prime*Q_S*M_prime/M/L+tau_H0;
    Q_z=L_prime*Q_W*Q_G;
    
    
    % initalize m_* with small values
    m_S=1e-3;
    m_G=1e-3;
    m_W=L_prime*m_S*M_prime/M/L+tau_H0;
    m_Z=L_prime*(m_W*m_G);
    
    %noise power
    SNR=snr(i);
    tau_N=10^(SNR/10*-1);
    
    for iii=1:maxite
        
        % \tilde m_Z in eq. (37a)
        mhat1=1/(tau_N+K*(Q_z-m_Z)*tau_X);
        m_zq=T*tau_X*mhat1;
        
        % \tilde m_W and \tilde m_G in eqs. (37b)--(37c)
        mhat2=1/(1/m_zq+Q_z-L_prime*m_W*m_G);
        m_wq=K*mhat2*m_G;
        m_gq=M*mhat2*m_W;
        
        % compute (37g) via numerical intergral
        pin=@(x,y) 1./(1+(1-lambda_G)/lambda_G*(m_gq+1)*exp(-(x.^2+y.^2)/(1+1/m_gq)));
        Uhat1=@(x,y) (x.^2+y.^2)/(m_gq+1)/(1+1/m_gq).*(pin(x,y).^2);
        pin=@(x,y) (1-lambda_G)/lambda_G*(m_gq+1)*exp(-(x.^2+y.^2)*m_gq);
        Uhat2=@(x,y) (x.^2+y.^2)/(1+1/m_gq)./((1+pin(x,y)).^2);
        term3=@(x,y) Uhat1(x,y).*exp(-(x.^2+y.^2))/(pi);
        term4=@(x,y) Uhat2(x,y).*exp(-(x.^2+y.^2))/(pi);
        m_G=(1-lambda_G)*integral2(term3,-int_region,int_region,-int_region,int_region)...
            +lambda_G*integral2(term4,-int_region,int_region,-int_region,int_region);
        
        
        % \tilde m_S in eq. (37d)
        mhat3=1/(1/m_wq+L_prime*(Q_S-m_S)*M_prime/M/L);
        m_sq=L_prime/L*mhat3;
        
        % compute (37h) via numerical intergral
        pin=@(x,y) 1./(1+(1-lambda_S)/lambda_S*(m_sq+1)*exp(-(x.^2+y.^2)/(1+1/m_sq)));
        Uhat1=@(x,y) (x.^2+y.^2)/(m_sq+1)/(1+1/m_sq).*(pin(x,y).^2);
        pin=@(x,y) (1-lambda_S)/lambda_S*(m_sq+1)*exp(-(x.^2+y.^2)*m_sq);
        Uhat2=@(x,y) (x.^2+y.^2)/(1+1/m_sq)./((1+pin(x,y)).^2);
        term3=@(x,y) Uhat1(x,y).*exp(-(x.^2+y.^2))/(pi);
        term4=@(x,y) Uhat2(x,y).*exp(-(x.^2+y.^2))/(pi);
        m_S=(1-lambda_S)*integral2(term3,-int_region,int_region,-int_region,int_region)+lambda_S*integral2(term4,-int_region,int_region,-int_region,int_region);
        
        
        
        % eq. (37f)
        v_w=(Q_W-tau_H0-L_prime*m_S*M_prime/M/L)/...
            (m_wq*(Q_W-tau_H0-L_prime*m_S*M_prime/M/L)+1);
        m_W=Q_W-v_w;
        % eq. (37e)
        v_z=(Q_z-L_prime*m_W*m_G)/(m_zq*(Q_z-L_prime*m_W*m_G)+1);
        m_Z=Q_z-v_z;
        
        
        % compute MSEs in dB by eqs. (37i)--(37j)
        eG=10*log10((Q_G-m_G));
        eS=10*log10((Q_S-m_S));
        
        %early stop
        if abs(MSE_S(i)-eS)<=threshold && abs(MSE_G(i)-eG)<=threshold
            MSE_S(i)=eS;
            MSE_G(i)=eG;
            fprintf('i: %d, v_z: %f m_zq:%f v_w: %f m_f: %f Er %f, m_x: %f Eu %f\n'...
                ,iii,10*log10(v_z),m_zq,10*log10(v_w),m_S,MSE_S(i),m_G,MSE_G(i));
            
            break;
        end
        if eS<= -90 && eG<=-90
            MSE_S(i)=eS;
            MSE_G(i)=eG;
            fprintf('i: %d, v_z: %f m_zq:%f v_w: %f m_f: %f Er %f, m_x: %f Eu %f\n'...
                ,iii,10*log10(v_z),m_zq,10*log10(v_w),m_S,MSE_S(i),m_G,MSE_G(i));
            break;
        end
        
        % store MSEs
        MSE_S(i)=eS;
        MSE_G(i)=eG;
    end
end


end
