function hplot = dibujar_criterio_dano1(ce,nu,q,tipo_linea,MDtype,n)
%*************************************************************************************
%*                 PLOT DAMAGE SURFACE CRITERIUM: ISOTROPIC MODEL                             %*
%*                                                                                  %*
%*      function [ce] = tensor_elastico (Eprop, ntype)                    %*
%*                                                                                  %*
%*      INPUTS                                                       %*
%*                                                                                  %*
%*                    Eprop(4)    vector de propiedades de material                 %*
%*                                      Eprop(1)=  E------>modulo de Young          %*
%*                                      Eprop(2)=  nu----->modulo de Poisson        %*
%*                                      Eprop(3)=  H----->modulo de Softening/hard. %*
%*                                      Eprop(4)=sigma_u----->tensi�n �ltima        %*
%*                     ntype                                 %*
%*                                 ntype=1  plane stress                            %*
%*                                 ntype=2  plane strain                            %*
%*                                 ntype=3  3D                                      %*
%*                     ce(4,4)     Constitutive elastic tensor  (PLANE S.       )    %*
%*                     ce(6,6)                                  ( 3D)                %*
%*************************************************************************************


%*************************************************************************************
%*        Inverse ce                                                                %*
ce_inv=inv(ce);
c11=ce_inv(1,1);
c22=ce_inv(2,2);
c12=ce_inv(1,2);
c21=c12;
c14=ce_inv(1,4);
c24=ce_inv(2,4);
%**************************************************************************************







%**************************************************************************************
% POLAR COORDINATES
if MDtype==1
    tetha=[0:0.01:2*pi];
    %**************************************************************************************
    %* RADIUS
    D=size(tetha);                       %*  Range
    m1=cos(tetha);                       %*
    m2=sin(tetha);                       %*
    Contador=D(1,2);                     %*
    
    
    radio = zeros(1,Contador) ;
    s1    = zeros(1,Contador) ;
    s2    = zeros(1,Contador) ;
    
    for i=1:Contador
        radio(i)= q/sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[m1(i) m2(i) 0 ...
            nu*(m1(i)+m2(i))]');
        
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);  
        
    end
    hplot =plot(s1,s2,tipo_linea);
    
    
elseif MDtype==2
 %   Comment/delete lines below once you have implemented this case
  %  *******************************************************
   tetha=[-pi/2*0.8:0.01:pi*0.9];
    %**************************************************************************************
    %* RADIUS
    D=size(tetha);                       %*  Range
    m1=cos(tetha);                       %*
    m2=sin(tetha);                       %*
    Contador=D(1,2);                     %*
    
    
    radio = zeros(1,Contador) ;
    s1    = zeros(1,Contador) ;
    s2    = zeros(1,Contador) ;
    
    for i=1:Contador
        
        radio(i)= q/sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[m1(i) m2(i) 0 ...
            nu*(m1(i)+m2(i))]');
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);
        
        if tetha(i)<0
        radio(i)= q/sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[m1(i) 0 0 ...
            nu*(m1(i)+0)]'); %m2(i)=0
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);
        end
        if tetha(i)>pi/2
        radio(i)= q/sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[0 m2(i) 0 ...
            nu*(0+m2(i))]'); %m1(i)=0
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);
        end

    end
   
    hplot =plot(s1,s2,tipo_linea);
     
  
  
    
elseif MDtype==3
    % Comment/delete lines below once you have implemented this case
    % *******************************************************
    tetha=[0:0.001:2*pi];
    %**************************************************************************************
    %* RADIUS
    D=size(tetha);                       %*  Range
    m1=cos(tetha);                       %*
    m2=sin(tetha);                       %*
    Contador=D(1,2);                     %*
    
    
    radio = zeros(1,Contador) ;
    s1    = zeros(1,Contador) ;
    s2    = zeros(1,Contador) ;

    for i=1:Contador
        
        theta2 = ((m1(i)>0)*m1(i) + (m2(i)>0)*m2(i)  + (nu*(m1(i)+m2(i)))*((nu*(m1(i)+m2(i)))>0)    )/((abs(m1(i))) + (abs(m2(i))) + (abs(nu*(m1(i)+m2(i)))));
        
        radio(i)= q/((theta2+(1-theta2)/n)*sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[m1(i) m2(i) 0 ...
            nu*(m1(i)+m2(i))]'));
        
        
        
        
        
        
%         m3(i)=nu*(m1(i)+m2(i));
%         m1_mc(i)=m1(i)*(m1(i)>0);
%         m2_mc(i)=m2(i)*(m2(i)>0);
%         m3_mc(i)=nu*(m1_mc(i)+m2_mc(i))*((nu*(m1_mc(i)+m2_mc(i))>0));
%         tetha2(i)=(m1_mc(i)+m2_mc(i)+m3_mc(i))/(m1(i)+m1(i)+m3(i));
%         radio(i)= (tetha2(i)+(1-tetha2(i))/n)*q/sqrt([m1(i) m2(i) 0 nu*(m1(i)+m2(i))]*ce_inv*[m1(i) m2(i) 0 ...
%             nu*(m1(i)+m2(i))]');
        
        s1(i)=radio(i)*m1(i);
        s2(i)=radio(i)*m2(i);  
        
    end
    hplot =plot(s1,s2,tipo_linea);
    
    
    
    
    
end
%**************************************************************************************



%**************************************************************************************
return



