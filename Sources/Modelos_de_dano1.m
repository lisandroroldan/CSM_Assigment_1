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
if (MDtype==1)      %* Symmetric
rtrial= sqrt(eps_n1*ce*eps_n1')                       ;

elseif (MDtype==2)  %* Only tension 

    
elseif (MDtype==3)  %*Non-symmetric

end
%**************************************************************************************
return