function [sigma_n1,hvar_n1,aux_var] = rmap_dano1 (eps_n1,hvar_n,Eprop,ce,MDtype,n)

%**************************************************************************************
%*                                         *
%*           Integration Algorithm for a isotropic damage model
%*
%*                                                                                    *
%*            [sigma_n1,hvar_n1,aux_var] = rmap_dano1 (eps_n1,hvar_n,Eprop,ce)        *
%*                                                                                    *
%* INPUTS              eps_n1(4)   strain (almansi)    step n+1                       *
%*                                 vector R4    (exx eyy exy ezz)                     *
%*                     hvar_n(6)   internal variables , step n                        *
%*                                 hvar_n(1:4) (empty)                          *
%*                                 hvar_n(5) = r  ; hvar_n(6)=q                       *
%*                     Eprop(:)    Material parameters                                *
%*
%*                     ce(4,4)     Constitutive elastic tensor                        *
%*                                                                                    *
%* OUTPUTS:            sigma_n1(4) Cauchy stress  , step n+1                          *
%*                     hvar_n(6)   Internal variables , step n+1                           *
%*                     aux_var(3)  Auxiliar variables for computing const. tangent tensor  *
%***************************************************************************************


hvar_n1 = hvar_n;
r_n     = hvar_n(5);
q_n     = hvar_n(6);
E       = Eprop(1);
nu      = Eprop(2);
H       = Eprop(3);
sigma_u = Eprop(4);
hard_type = Eprop(5) ;
%*************************************************************************************


%*************************************************************************************
%*       initializing                                                %*
 r0 = sigma_u/sqrt(E);
 zero_q=1.d-6*r0;
% if(r_n<=0.d0)
%     r_n=r0;
%     q_n=r0;
% end
%*************************************************************************************


%*************************************************************************************
%*       Damage surface                                                              %*
[rtrial] = Modelos_de_dano1 (MDtype,ce,eps_n1,n);
%*************************************************************************************


%*************************************************************************************
%*   Ver el Estado de Carga                                                           %*
%*   --------->    fload=0 : elastic unload                                           %*
%*   --------->    fload=1 : damage (compute algorithmic constitutive tensor)         %*
fload=0;

if(rtrial > r_n)
    %*   Loading

    fload=1;
    delta_r=rtrial-r_n;
    r_n1= rtrial  ;
    if hard_type == 0
        %  Linear
        q_n1= q_n+ H*delta_r;
    else
        % Comment/delete lines below once you have implemented this case
        % *******************************************************
        menu({'Hardening/Softening exponential law has not been implemented yet. '; ...
            'Modify file "rmap_dano1" ' ; ...
            'to include this option'},  ...
            'STOP');
        error('OPTION NOT AVAILABLE')
    end

    if(q_n1<zero_q)
        q_n1=zero_q;
    end


else

    %*     Elastic load/unload
    fload=0;
    r_n1= r_n  ;
    q_n1= q_n  ;


end
% Damage variable
% ---------------
dano_n1   = 1.d0-(q_n1/r_n1);
%  Computing stress
%  ****************
sigma_n1  =(1.d0-dano_n1)*ce*eps_n1';
%hold on 
%plot(sigma_n1(1),sigma_n1(2),'bx')

%*************************************************************************************


%*************************************************************************************
%* Updating historic variables                                            %*
%  hvar_n1(1:4)  = eps_n1p;
hvar_n1(5)= r_n1 ;
hvar_n1(6)= q_n1 ;
%*************************************************************************************




%*************************************************************************************
%* Auxiliar variables                                                               %*
aux_var(1) = fload;
aux_var(2) = q_n1/r_n1;
%*aux_var(3) = (q_n1-H*r_n1)/r_n1^3;
%*************************************************************************************
 









