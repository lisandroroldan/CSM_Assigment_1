function [ce] = tensor_elastico1 (Eprop, ntype)
%*************************************************************************************
%*       Elastic constitutive tensor                                                %*
%************************************************************************************* 


%*************************************************************************************
%
%*                   G -------->  Shear modulus                                     %*                          
%*                   K -------->  Bulk modulus                                      %*
G=Eprop(1)/(2*(1+Eprop(2)));
K=Eprop(1)/(3*(1-2*Eprop(2)));
%*************************************************************************************


%*************************************************************************************
if(ntype==1)                            % Plane stress      
elseif(ntype==2)                        % Plane strain 
        ce    = zeros(4,4);             % Init.
        C1=K+(4.0D0/3.0D0)*G;           
        C2=K-(2.0D0/3.0D0)*G;         
        ce(1,1)=C1;                                                                 
        ce(2,2)=C1;
        ce(4,4)=C1;
        ce(1,2)=C2;
        ce(1,4)=C2;
        ce(2,4)=C2;
        ce(2,1)=C2;
        ce(4,1)=C2;
        ce(4,2)=C2;
        ce(3,3)=G;
elseif(ntype==4)                        % Tres Dimensiones
end
%*************************************************************************************

return