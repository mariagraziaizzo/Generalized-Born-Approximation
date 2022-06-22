clear all
close all
format long e
%
n0=2.292; % longitudinal phase velocity INPUT PARAMETER
rlt=0.54; % ratio of transverse to longitudinal phase velocity INPUT PARAMETER
n0t=rlt*n0; % transverse phase velocity
n0T=n0t;
l=0.143; %correlation length (a) INPUT PARAMETER
A=1.825;% intensity of fluctuations of elastic constant (including normalizing factors!) INPUT PARAMETER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Art=6; % strenght of coupling with intermolecular vibrational modes with transverse dynamics INPUT PARAMETER Art=0 no coupling
Arl=20; % strenght of coupling with intermolecular vibrational modes with longitudinal dynamics INPUT PARAMETER Art=0 no coupling
wr=7.2; % characteristic frequency of IVM INPUT PARAMETER
Gr=7.5; % damping of IVM INPUT PARAMETER
% 
Art_1=4; % as above, this is the second IVM
Arl_1=1; 
wr_1=11; 
Gr_1=6.5;
%
n_0T=(n0t^2-A*(2/5)*(1/n0t^2+(2/3)/n0^2))^0.5
n_0t=n_0T;
n_0L=(n0^2-A*(4/5)*(1/n0^2+(2/3)/n0t^2))^0.5
A1=A*(n_0T^4/n0t^4);
A1_L=A*(n_0T^2/n0t^2);
%
kD=17*0.62;
Qm=17.9;
Q=(5:0.1:30);
dQm=3.2;
s=100;
k=(0.1:0.3:16);
step=0.002;
first=0.01;
last=95;
step2=0.0005; % used to cut the tails
x_data=(first:step:last);
w=zeros(length(x_data),length(k));
for j=1:length(k);
w(:,j)=x_data;
end
wx=(-s:step:s);
w_L=zeros(length(wx),length(k));
for j=1:length(k);
    w_L(:,j)=wx;
end
L=zeros(length(wx),length(k));
LT=zeros(length(wx),length(k));
xConvol=(first-s:step:last+s);
w_Convol=zeros(length(xConvol),length(k));
for j=1:length(k);
w_Convol(:,j)=xConvol;
end
e=0.00000000001;
S0L=zeros(length(x_data),length(k));
C0L=zeros(length(x_data),length(k));
%
%%%%%%%%%%%%%%%%%%% LONGITUDINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(k);
S0L(:,j)=x_data./abs(x_data).*imag(1./(-((x_data+i*e).^2)./(n0^2)+k(j).^2));
C0L(:,j)=(x_data./abs(x_data).*imag(1./(-((x_data+i*e).^2)./(n0^2)+k(j).^2))).*(x_data.^2)./(k(j).^2);
S0L(:,j)=S0L(:,j)./max(S0L(:,j));
C0L(:,j)=C0L(:,j)./max(C0L(:,j));
end
% %%%%%%%%%%%%%%%% TRANSVERSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S0T=zeros(length(x_data),length(k));
C0T=zeros(length(x_data),length(k));
for j=1:length(k);
S0T(:,j)=x_data./abs(x_data).*imag(1./(-((x_data+i*e).^2)./(n0t^2)+k(j).^2));
C0T(:,j)=(x_data./abs(x_data).*imag(1./(-((x_data+i*e).^2)./(n0t^2)+k(j).^2))).*(x_data.^2./k(j).^2);
S0T(:,j)=S0T(:,j)./max(S0T(:,j));
C0T(:,j)=C0T(:,j)./max(C0T(:,j));
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0ll=zeros(length(x_data),length(k));
A1ll=zeros(length(x_data),length(k));
A2ll=zeros(length(x_data),length(k));
C1ll=zeros(length(x_data),length(k));
C2ll=zeros(length(x_data),length(k));
D1ll=zeros(length(x_data),length(k));
D2ll=zeros(length(x_data),length(k));
D1tt=zeros(length(x_data),length(k));
D2tt=zeros(length(x_data),length(k));
D1lt=zeros(length(x_data),length(k));
D2lt=zeros(length(x_data),length(k));
D1tl=zeros(length(x_data),length(k));
D2tl=zeros(length(x_data),length(k));
F1ll=zeros(length(x_data),length(k));
F2ll=zeros(length(x_data),length(k));
F1tt=zeros(length(x_data),length(k));
F2tt=zeros(length(x_data),length(k));
F1lt=zeros(length(x_data),length(k));
F2lt=zeros(length(x_data),length(k));
F1tl=zeros(length(x_data),length(k));
F2tl=zeros(length(x_data),length(k));
E0lt=zeros(length(x_data),length(k));
A1lt=zeros(length(x_data),length(k));
A2lt=zeros(length(x_data),length(k));
C1lt=zeros(length(x_data),length(k));
C2lt=zeros(length(x_data),length(k));
E0tt=zeros(length(x_data),length(k));
A1tt=zeros(length(x_data),length(k));
A2tt=zeros(length(x_data),length(k));
C1tt=zeros(length(x_data),length(k));
C2tt=zeros(length(x_data),length(k));
E0tl=zeros(length(x_data),length(k));
A1tl=zeros(length(x_data),length(k));
A2tl=zeros(length(x_data),length(k));
C1tl=zeros(length(x_data),length(k));
C2tl=zeros(length(x_data),length(k));
E0=zeros(length(x_data),length(k));
E0T=zeros(length(x_data),length(k));
E0_mu=zeros(length(x_data),length(k));
E0T_mu=zeros(length(x_data),length(k));
p0L1=zeros(length(x_data),length(k));
p0T1=zeros(length(x_data),length(k));
%
step=0.01;
t=(-1:step:1);
g=(-1:step:1);
%
%%%%%%%%%%%%%%%%%%%%%%%%% LONGITUDINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(k);
    for d=1:length(t);
    A2ll(:,j)=4.*0.5.*(t(d).^4).*((w(:,j)./n_0L).^3).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0L).^2+2.*k(j).*t(d).*w(:,j)./n_0L).^2);
    A1ll(:,j)=A1ll(:,j)+A2ll(:,j).*step;
    l1=(l^(-2)+k(j).^2.*(1-t(d).^2)).^0.5;
    C2ll(:,j)=-4.*((1+aFS.*abs(t(d))).*(t(d).^4).*((k(j)*t(d)+i*l1).^3)./(2.*l1.^2)).*(1./((w(:,j)./n_0L).^2-(k(j).*t(d)+i*l1).^2)).*(2-((k(j).*t(d)+i*l1)./(2*i*l1))-(((k(j).*t(d)+i*l1).^2)./((k(j).*t(d)+i*l1).^2-((w(:,j)./n_0L).^2))));
    C1ll(:,j)=C1ll(:,j)+C2ll(:,j).*step;
    A2lt(:,j)=4.*0.5.*(t(d).^2-t(d).^4).*((w(:,j)./n_0t).^3).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0t).^2+2.*k(j).*t(d).*w(:,j)./n_0t).^2);
    A1lt(:,j)=A1lt(:,j)+A2lt(:,j).*step;
    C2lt(:,j)=-4.*((1+aFS.*abs(t(d))).*(t(d).^2-t(d).^4).*((k(j)*t(d)+i*l1).^3)./(2.*l1.^2)).*(1./((w(:,j)./n_0t).^2-(k(j).*t(d)+i*l1).^2)).*(2-((k(j).*t(d)+i*l1)./(2*i*l1))-(((k(j).*t(d)+i*l1).^2)./((k(j).*t(d)+i*l1).^2-((w(:,j)./n_0t).^2))));
    C1lt(:,j)=C1lt(:,j)+C2lt(:,j).*step;
    D2ll(:,j)=-4.*(t(d).^4).*((w(:,j)./n_0L).^3./4).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0L).^2+2.*k(j).*(w(:,j)./n_0L).*t(d)).^2).*(5-4.*(w(:,j)./n_0L).*(w(:,j)./n_0L+k(j).*t(d))./(l^(-2)+k(j).^2+(w(:,j)./n_0L).^2+2.*k(j).*t(d).*(w(:,j)./n_0L)));
    D1ll(:,j)=D1ll(:,j)+D2ll(:,j).*step; %(Res(-p0))
    F2ll(:,j)=-(t(d).^4).*(1./l1.^2).*((k(j).*t(d)+i.*l1).^5./((w(:,j)./n_0L).^2-(k(j).*t(d)+i.*l1).^2).^2).*(5+4.*(k(j).*t(d)+i.*l1).^2./((w(:,j)./n_0L).^2-(k(j).*t(d)+i.*l1).^2)-k(j).*t(d)./(i.*l1));
    F1ll(:,j)=F1ll(:,j)+F2ll(:,j).*step; % Res(kx+ia(k,x))
    D2lt(:,j)=-4.*(t(d).^2-t(d).^4).*((w(:,j)./n_0t).^3./4).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0t).^2+2.*k(j).*(w(:,j)./n_0t).*t(d)).^2).*(5-4.*(w(:,j)./n_0t).*(w(:,j)./n_0t+k(j).*t(d))./(l^(-2)+k(j).^2+(w(:,j)./n_0t).^2+2.*k(j).*t(d).*(w(:,j)./n_0t)));
    D1lt(:,j)=D1lt(:,j)+D2lt(:,j).*step; %(Res(-p0))
    F2lt(:,j)=-(t(d).^2-t(d).^4).*(1./l1.^2).*((k(j).*t(d)+i.*l1).^5./((w(:,j)./n_0t).^2-(k(j).*t(d)+i.*l1).^2).^2).*(5+4.*(k(j).*t(d)+i.*l1).^2./((w(:,j)./n_0t).^2-(k(j).*t(d)+i.*l1).^2)-k(j).*t(d)./(i.*l1));
    F1lt(:,j)=F1lt(:,j)+F2lt(:,j).*step; % Res(kx+ia(k,x))
    end
    E0ll(:,j)=i*(A/(n_0L^2)).*((k(j)^2)./(1)).*(1/l).*(2.*A1ll(:,j)+2.*C1ll(:,j))+i.*((A^2/(n_0L^4))).*(((k(j)).^2)./(1)).*((1/l).^2).*(2.*D1ll(:,j)+2.*F1ll(:,j)).*((1/n0^2).*(4/5).*(1./((w(:,j)./n0).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0).^2+l^(-2))+i.*2.*(w(:,j)./n0).^3+l^(-3))+(1/n0t^2).*(8/15).*(1./((w(:,j)./n0t).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0t).^2+l^(-2))+i.*2.*(w(:,j)./n0t).^3+l^(-3)));
    E0lt(:,j)=i*(A/(n_0t^2)).*((k(j)^2)./(1)).*(1/l).*(2.*A1lt(:,j)+2.*C1lt(:,j))+i.*((A^2/(n_0t^4))).*(((k(j)).^2)./(1)).*((1/l).^2).*(2.*D1lt(:,j)+2.*F1lt(:,j)).*0.5.*((1/n0t^2).*(4/5).*(1./((w(:,j)./n0t).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0t).^2+l^(-2))+i.*2.*(w(:,j)./n0t).^3+l^(-3))+(1/n0^2).*(8/15).*(1./((w(:,j)./n0).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0).^2+l^(-2))+i.*2.*(w(:,j)./n0).^3+l^(-3)));
    E0_mu(:,j)=E0ll(:,j)+E0lt(:,j);%E0l=E0ll+E0lt 
    E0(:,j)=E0_mu(:,j)+(Arl.*k(j).^2)./(w(:,j).^2-wr^2-i*Gr.*w(:,j))+(Arl_1.*k(j).^2)./(w(:,j).^2-wr_1^2-i*Gr_1.*w(:,j));
end
I=exp(-(Q-Qm).^2./(dQm^2));
I=I./trapz(I);
%
SL=zeros(length(x_data),length(k));
SL_bend=zeros(length(x_data),length(k));
CL=zeros(length(x_data),length(k));
CL_bend=zeros(length(x_data),length(k));
SL_bend_d=zeros(length(x_data),length(k));
CL_bend_d=zeros(length(x_data),length(k));
SL_bend_dh=zeros(length(x_data),length(k),length(Q));
SL_d=zeros(length(x_data),length(k));
CL_d=zeros(length(x_data),length(k));
SL_dh=zeros(length(x_data),length(k),length(Q));
for j=1:length(k);
r=k(j);    
SL(:,j)=(k(j).^2).*(abs(imag(E0(:,j)))./w(:,j))./(((w(:,j).^2)-(n0.*k(j))^2-real(E0(:,j))).^2+imag(E0(:,j)).^2);
SL_bend(:,j)=(k(j).^2).*(abs(imag(E0(:,j)))./(w(:,j)./((1/r).*(Qm./pi).*abs(sin(r.*pi./Qm)))))./((((w(:,j)./((1/r).*(Qm./pi).*abs(sin(r.*pi./Qm)))).^2)-(n0.*k(j))^2-real(E0(:,j))).^2+imag(E0(:,j)).^2); % dynamic structure factor including the bending

for h=1:length(Q)
SL_bend_dh(:,j,h)=I(h).*(k(j).^2).*(abs(imag(E0(:,j)))./(w(:,j)./(((Q(h)./pi).*abs(sin(r.*pi./Q(h))))./((Qm./pi).*abs(sin(r.*pi./Qm))))))./((((w(:,j)./(((Q(h)./pi).*abs(sin(r.*pi./Q(h))))./((Qm./pi).*abs(sin(r.*pi./Qm))))).^2)-(n0.*k(j))^2-real(E0(:,j))).^2+imag(E0(:,j)).^2); %.*(1/r).*(Q./pi).*abs(sin(r.*pi./Q))
SL_dh(:,j,h)=I(h).*(k(j).^2).*(abs(imag(E0(:,j)))./w(:,j))./(((w(:,j).^2)-(n0.*(((Q(h)/pi).*abs(sin(k(j).*pi./Q(h))))./((Qm/pi).*abs(sin(k(j).*pi./Qm)))).*k(j)).^2-real(E0(:,j))).^2+imag(E0(:,j)).^2);
end
SL_bend_d(:,j)=sum(SL_bend_dh(:,j,:),3);
SL_d(:,j)=sum(SL_dh(:,j,:),3);
SL(:,j)=SL(:,j)./max(SL(:,j));
SL_bend(:,j)=SL_bend(:,j)./max(SL_bend(:,j));
SL_bend_d(:,j)=SL_bend_d(:,j)./max(SL_bend_d(:,j));
SL_d(:,j)=SL_d(:,j)./max(SL_d(:,j));
CL(:,j)=SL(:,j).*(w(:,j).^2./k(j)^2);
CL(:,j)=CL(:,j)./max(CL(:,j));
CL_bend(:,j)=SL_bend(:,j).*(w(:,j).^2./k(j)^2);
CL_bend(:,j)=CL_bend(:,j)./max(CL_bend(:,j));
CL_bend_d(:,j)=SL_bend_d(:,j).*(w(:,j).^2./k(j)^2);
CL_bend_d(:,j)=CL_bend_d(:,j)./max(CL_bend_d(:,j));
CL_d(:,j)=SL_d(:,j).*(w(:,j).^2./k(j)^2);
CL_d(:,j)=CL_d(:,j)./max(CL_d(:,j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dw_p=(n_0L*(Qm+dQm)/pi).*abs(sin(k.*pi./(Qm+dQm)));
dw_n=(n_0L*(Qm-dQm)/pi).*abs(sin(k.*pi./(Qm-dQm)));
dw_pT=(n_0T*(Qm+dQm)/pi).*abs(sin(k.*pi./(Qm+dQm)));
dw_nT=(n_0T*(Qm-dQm)/pi).*abs(sin(k.*pi./(Qm-dQm)));
%
figure(1)
plot(Q,I);
%
figure(2)
plot(k,dw_p,'* green');
hold on
plot(k,dw_n,'* yellow');
%
figure(3)
plot(k,dw_pT,'* green');
hold on
plot(k,dw_nT,'* yellow');
%
figure(4)
plot(k,dw_p./((Qm/pi).*abs(sin(k.*pi./Qm))),'green *');
hold on
plot(k,dw_n./((Qm/pi).*abs(sin(k.*pi./Qm))),'yellow *');
%
figure(5)
plot(k,dw_pT./((Qm/pi).*abs(sin(k.*pi./Qm))),'green *');
hold on
plot(k,dw_nT./((Qm/pi).*abs(sin(k.*pi./Qm))),'yellow *');
%
figure(6)
G=0.5.*abs(dw_p-dw_n); %HWHM
plot(k,G,'red');
%
figure(7)
plot(k,abs((dw_pT./((Qm/pi).*abs(sin(k.*pi./Qm))))-(dw_nT./((Qm/pi).*abs(sin(k.*pi./Qm))))));
hold on
GT=0.5.*abs(dw_pT-dw_nT); %HWHM
plot(k,GT,'red');
%
for j=1:length(k)
    L(:,j)=exp(-w_L(:,j).^2./((2/1.4)*(G(j)^2)));
    I_L=trapz(L(:,j));
    L(:,j)=L(:,j)./I_L;
    L(:,j)=L(:,j).*step;
    LT(:,j)=exp(-w_L(:,j).^2./(GT(j)^2));
    I_LT=trapz(LT(:,j));
    LT(:,j)=LT(:,j)./I_LT;
    LT(:,j)=LT(:,j).*step;
end
%
figure(8)
for j=1:length(k)
    plot(wx,L(:,j)./max(L(:,j)));
    hold on
end
%
figure(9)
for j=1:length(k)
    plot(wx,LT(:,j)./max(LT(:,j)));
    hold on
end
%
SL_c=zeros(length(x_data),length(k));
CL_c=zeros(length(x_data),length(k));
Convol_L=zeros(length(x_data),length(k));
SL_bend_c=zeros(length(x_data),length(k));
CL_bend_c=zeros(length(x_data),length(k));
Convol_L_bend=zeros(length(x_data),length(k));
Convol_CL_bend=zeros(length(x_data),length(k));
%
for j=1:length(k);
    r=k(j);
Lj=L(:,j);
SLj=SL(:,j);
SL_bendj=SL_bend(:,j);
CL_bendj=CL_bend(:,j);
Convol_Lj=conv(Lj,SLj);
Convol_L_bendj=conv(Lj,SL_bendj);
Convol_CL_bendj=conv(Lj,CL_bendj);
SL_c(:,j)= interp1(w_Convol(:,j),Convol_Lj,w(:,j));
SL_bend_c(:,j)= interp1(w_Convol(:,j),Convol_L_bendj,w(:,j));
CL_bend_c(:,j)= interp1(w_Convol(:,j),Convol_CL_bendj,w(:,j));
CL_c(:,j)=SL_c(:,j).*(w(:,j).^2)./(k(j).^2);
SL_c(:,j)=SL_c(:,j)./max(SL_c(:,j));
CL_c(:,j)=CL_c(:,j)./max(CL_c(:,j));
SL_bend_c(:,j)=SL_bend_c(:,j)./max(SL_bend_c(:,j));
CL_bend_c(:,j)=CL_bend_c(:,j)./max(CL_bend_c(:,j));
end
% 
%%%%%%%%%%%%%%%%%%%%%% TRANSVERSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(k);
    for d=1:length(t);
      A2tt(:,j)=0.5.*(1-3.*t(d).^2+4.*t(d).^4).*((w(:,j)./n_0t).^3).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0t).^2+2.*k(j).*t(d).*w(:,j)./n_0t).^2);
      A1tt(:,j)=A1tt(:,j)+A2tt(:,j).*step;
      l1=(l^(-2)+k(j).^2.*(1-t(d).^2)).^0.5;
      C2tt(:,j)=-((1-3.*t(d).^2+4.*t(d).^4).*((k(j)*t(d)+i*l1).^3)./(2.*l1.^2)).*(1./((w(:,j)./n_0t).^2-(k(j).*t(d)+i*l1).^2)).*(2-((k(j).*t(d)+i*l1)./(2*i*l1))-(((k(j).*t(d)+i*l1).^2)./((k(j).*t(d)+i*l1).^2-((w(:,j)./n_0t).^2))));
      C1tt(:,j)=C1tt(:,j)+C2tt(:,j).*step;
      A2tl(:,j)=0.5.*4.*(t(d).^2-t(d).^4).*((w(:,j)./n_0L).^3).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0L).^2+2.*k(j).*t(d).*w(:,j)./n_0L).^2);
      A1tl(:,j)=A1tl(:,j)+A2tl(:,j).*step;
      l1=(l^(-2)+k(j).^2.*(1-t(d).^2)).^0.5;
      C2tl(:,j)=-4.*((t(d).^2-t(d).^4).*((k(j)*t(d)+i*l1).^3)./(2.*l1.^2)).*(1./((w(:,j)./n_0L).^2-(k(j).*t(d)+i*l1).^2)).*(2-((k(j).*t(d)+i*l1)./(2*i*l1))-(((k(j).*t(d)+i*l1).^2)./((k(j).*t(d)+i*l1).^2-((w(:,j)./n_0L).^2))));
      C1tl(:,j)=C1tl(:,j)+C2tl(:,j).*step;
      D2tt(:,j)=-(1-3.*t(d).^2+4.*t(d).^4).*((w(:,j)./n_0t).^3./4).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0t).^2+2.*k(j).*(w(:,j)./n_0t).*t(d)).^2).*(5-4.*(w(:,j)./n_0t).*(w(:,j)./n_0t+k(j).*t(d))./(l^(-2)+k(j).^2+(w(:,j)./n_0t).^2+2.*k(j).*t(d).*(w(:,j)./n_0t)));
      D1tt(:,j)=D1tt(:,j)+D2tt(:,j).*step; %(Res(-p0))
      F2tt(:,j)=-(1/4).*(1-3.*t(d).^2+4.*t(d).^4).*(1./l1.^2).*((k(j).*t(d)+i.*l1).^5./((w(:,j)./n_0t).^2-(k(j).*t(d)+i.*l1).^2).^2).*(5+4.*(k(j).*t(d)+i.*l1).^2./((w(:,j)./n_0t).^2-(k(j).*t(d)+i.*l1).^2)-k(j).*t(d)./(i.*l1));
      F1tt(:,j)=F1tt(:,j)+F2tt(:,j).*step; % Res(kx+ia(k,x))
      D2tl(:,j)=-4.*(t(d).^2-t(d).^4).*((w(:,j)./n_0L).^3./4).*(1./(l^(-2)+k(j).^2+(w(:,j)./n_0L).^2+2.*k(j).*(w(:,j)./n_0L).*t(d)).^2).*(5-4.*(w(:,j)./n_0L).*(w(:,j)./n_0L+k(j).*t(d))./(l^(-2)+k(j).^2+(w(:,j)./n_0L).^2+2.*k(j).*t(d).*(w(:,j)./n_0L)));
      D1tl(:,j)=D1tl(:,j)+D2tl(:,j).*step; %(Res(-p0))
      F2tl(:,j)=-(t(d).^2-t(d).^4).*(1./l1.^2).*((k(j).*t(d)+i.*l1).^5./((w(:,j)./n_0L).^2-(k(j).*t(d)+i.*l1).^2).^2).*(5+4.*(k(j).*t(d)+i.*l1).^2./((w(:,j)./n_0L).^2-(k(j).*t(d)+i.*l1).^2)-k(j).*t(d)./(i.*l1));
      F1tl(:,j)=F1tl(:,j)+F2tl(:,j).*step; % Res(kx+ia(k,x))
    end
      E0tt(:,j)=i*(A/n_0t^2).*((k(j).^2)./(1)).*(1/l).*(2.*A1tt(:,j)+2.*C1tt(:,j))+i.*((A^2/(n_0t^4))).*(((k(j)).^2)./(1)).*((1/l).^2).*(2.*D1tt(:,j)+2.*F1tt(:,j)).*0.5.*((1/n0t^2).*(4/5).*(1./((w(:,j)./n0t).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0t).^2+l^(-2))+i.*2.*(w(:,j)./n0t).^3+l^(-3))+(1/n0^2).*(8/15).*(1./((w(:,j)./n0).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0).^2+l^(-2))+i.*2.*(w(:,j)./n0).^3+l^(-3)));
      E0tl(:,j)=i*(A/n_0L^2).*((k(j).^2)./(1)).*(1/l).*(2.*A1tl(:,j)+2.*C1tl(:,j))+i.*((A^2/(n_0L^4))).*(((k(j)).^2)./(1)).*((1/l).^2).*(2.*D1tl(:,j)+2.*F1tl(:,j)).*((1/n0^2).*(4/5).*(1./((w(:,j)./n0).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0).^2+l^(-2))+i.*2.*(w(:,j)./n0).^3+l^(-3))+(1/n0t^2).*(8/15).*(1./((w(:,j)./n0t).^2+l^(-2)).^2).*(-l^(-1).*(3.*(w(:,j)./n0t).^2+l^(-2))+i.*2.*(w(:,j)./n0t).^3+l^(-3)));
      E0T_mu(:,j)=(0.5.*E0tl(:,j)+0.5*E0tt(:,j));
      E0T(:,j)=E0T_mu(:,j)+(Art.*k(j).^2)./(w(:,j).^2-wr^2-i*Gr.*w(:,j))+(Art_1.*k(j).^2)./(w(:,j).^2-wr_1^2-i*Gr_1.*w(:,j));
end
%
ST=zeros(length(x_data),length(k));
ST_bend=zeros(length(x_data),length(k));
CT=zeros(length(x_data),length(k));
CT_bend=zeros(length(x_data),length(k));
ST_bend_c=zeros(length(x_data),length(k));
CT_bend_c=zeros(length(x_data),length(k));
Convol_T_bend=zeros(length(x_data),length(k));
%
for j=1:length(k);
 r=k(j);
 ST(:,j)=(k(j).^2).*(abs(imag(E0T(:,j)))./w(:,j))./(((w(:,j).^2)-(n0t.*k(j))^2-real(E0T(:,j))).^2+imag(E0T(:,j)).^2);
 ST_bend(:,j)=(k(j).^2).*(abs(imag(E0T(:,j)))./(w(:,j)./((1/r).*(Qm./pi).*abs(sin(r.*pi./Qm)))))./((((w(:,j)./((1/r).*(Qm./pi).*abs(sin(r.*pi./Qm)))).^2)-(n0t.*k(j))^2-real(E0T(:,j))).^2+imag(E0T(:,j)).^2); %.*(1/r).*(Q./pi).*abs(sin(r.*pi./Q))
 ST(:,j)=ST(:,j)./max(ST(:,j));
 ST_bend(:,j)=ST_bend(:,j)./max(ST_bend(:,j));
 CT(:,j)=ST(:,j).*(w(:,j).^2./k(j)^2);
 CT_bend(:,j)=ST_bend(:,j).*((w(:,j)./((1/r).*(Qm./pi).*abs(sin(r.*pi./Qm)))).^2./k(j)^2);
 CT(:,j)=CT(:,j)./max(CT(:,j));
 CT_bend(:,j)=CT_bend(:,j)./max(CT_bend(:,j));
end
%
ST_c=zeros(length(x_data),length(k));
CT_c=zeros(length(x_data),length(k));
Convol_T=zeros(length(x_data),length(k));
%
for j=1:length(k);
r=k(j);
Lj=LT(:,j);
STj=ST(:,j);
ST_bendj=ST_bend(:,j);
Convol_T_bendj=conv(Lj,ST_bendj);
Convol_Tj=conv(Lj,STj);
ST_c(:,j)= interp1(w_Convol(:,j),Convol_Tj,w(:,j));
ST_bend_c(:,j)= interp1(w_Convol(:,j),Convol_T_bendj,w(:,j));
CT_c(:,j)=ST_c(:,j).*(w(:,j).^2)./(k(j).^2);
ST_c(:,j)=ST_c(:,j)./max(ST_c(:,j));
CT_c(:,j)=CT_c(:,j)./max(CT_c(:,j));
CT_bend_c(:,j)=ST_bend_c(:,j).*((w(:,j)./((1/k(j)).*(Qm./pi).*abs(sin(r.*pi./Qm)))).^2./k(j)^2);
ST_bend_c(:,j)=ST_bend_c(:,j)./max(ST_bend_c(:,j));
CT_bend_c(:,j)=CT_bend_c(:,j)./max(CT_bend_c(:,j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% SAVING
%
save('CL.mat','CL');
save('CL_c.mat','CL_c');
save('CL_bend.mat','CL_bend');
save('CL_bend_c.mat','CL_bend_c');
save('CL_bend_d.mat','CL_bend_d');
save('SL.mat','SL');
save('SL_c.mat','SL_c');
save('SL_bend.mat','SL_bend');
save('SL_bend_c.mat','SL_bend_c');
save('SL_bend_d.mat','SL_bend_d');
save('S0L.mat','S0L');
save('C0L.mat','C0L');
save('S0T.mat','S0T');
save('C0T.mat','C0T');
