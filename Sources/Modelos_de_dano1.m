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
sigma_v=eps_n1*ce; %stress tensor
for i=1:4
    if sigma_v(i)<0
        sigma_v(i)=0;
    end
end
rtrial= sqrt(sigma_v*eps_n1'); %return value


end
    
%elseif (MDtype==3)  %*Non-symmetric

%end
%**************************************************************************************
