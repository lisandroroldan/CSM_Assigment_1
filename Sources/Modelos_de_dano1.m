function [rtrial] = Modelos_de_dano1 (MDtype,ce,eps_n1,n)
%**************************************************************************************
%*          Defining damage criterion surface                                        %*
%*                                                                                   %*
%*
%*                          MDtype=  1      : SYMMETRIC                              %*
%*                          MDtype=  2      : ONLY TENSION                           %*
%*                          MDtype=  3      : NON-SYMMETRIC                          %*
%*                                                                                   %*
%*                                                                                   %*
%* OUTPUT:                                                                           %*
%*                          rtrial                                                   %*               
%**************************************************************************************



%**************************************************************************************
switch MDtype
    
case 1 %* Symmetric
     
rtrial= sqrt(eps_n1*ce*eps_n1');

case 2  %* Only tension 
sigma_v=eps_n1*ce %stress tensor (x,y)
sigma_m=[sigma_v(1) sigma_v(4); sigma_v(4) sigma_v(2)] %voig notation --> tensor
[Q,sigma_p]=eigs(sigma_m) %find the diagonal components (1,2)
%McAuley bracket
for i=1:2
    for j=1:2
        if sigma_p(i,j)<0
            sigma_p(i,j)=0;
        end
    end
end
sigma_bar_m=Q*sigma_p*Q'; %going back to cartesian (x,y)
sigma_bar_v=[sigma_bar_m(1,1) sigma_bar_m(2,1) 0 sigma_bar_m(2,2)]; %tensor --> voig notation
rtrial= sqrt(sigma_bar_v*eps_n1') %return value


end
    
%elseif (MDtype==3)  %*Non-symmetric

%end
%**************************************************************************************
