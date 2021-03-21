%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
% Nguoi thuc hien: Luu Truong Khanh                                       %  
% Ngay kiem tra: 3-7-2002                                                 %  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>CHUONG TRINH<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%
% chuong trinh tinh bai toan dong luc hoc theo thoi gian t tam chu nhat   %   
% tren nen dan hoi bien tu do giai phuong trinh vi phan dao ham rieng:    %  
%                                                                         %  
%               M * x'' + K * x' +  C * x = R (1.1)                       %  
%                                                                         %  
%   bang 3 phuong phap so trong co hoc ket cau:                           %  
%   - Phuong phap phan tich dong luc hoc Newmark                          %  
%       + Voi hai kieu tai trong dieu hoa = 1                             %  
%       + Voi hai kieu tai trong song nen = 1                             %  
%   - Phuong phap phan tich tinh luc hoc khi chiu tai trong tinh          %  
%       + Dau vao:  iel, nel, apt, bpt, iel, nel, x, y, ff, kkt, mmt      %  
%                                                                         %  
%       + Dau ra:   w, noiluctg, noilucnut_tinh, dispt                    %  
%                                                                         %  
%       + Chuong trinh con:                                               %  
%         [kkt,mmt,ff] = khudkbiendaodongt(kkt,mmt,ff,bcdof,bcval)        %  
%                                                                         %      
% Tinh gia tri rieng va vec to rieng cua ma tran do cung va khoi luong    %  
%                                                                         %  
% tinh chuyen vi voi tai trong song nen thay doi tu 15000 len 1000        %  
%=========================================================================%      
% Giai doan 1                                                             %  
% 1 : nt1 - giai doan tai trong tang tuyen tinh                           %  
%=========================================================================%
% Giai doan 2                                                             %  
% 1 : nt1 - giai doan tai trong tang deu                                  %  
% 1 : nt2 - giai doan tai trong giam deu                                  %  
%=========================================================================%
% Giai doan 3                                                             %  
% 1 : nt1 - giai doan tai trong tang deu                                  %  
% 1 : nt2 - giai doan tai trong giam deu                                  %  
% 1 : nt3 - giai doan tai trong di ngang                                  %  
%=========================================================================%
% Tai trong song nen phan bo deu cuongdo q theo qui luat:                 %  
% 0 =< t =< teta1: f(t)=t/teta1 -  taitrong tang tuyen tinh               %  
% teta1=< t =< teta2: f(t)=1-t/teta2 - taitrong giam tuyen tinh           %  
% t > (teta1+teta2): f(t)=0 - taitrong bang khong                         %  
% Voi taitrong dieu hoa phanbo deu: q(x,y,t)=q*sin(omega*t)               %  
%=========================================================================%
% theo chu vi voi kieu phantu chu nhat 4 nut, 12 chuyen vi nut            %  
% thietlap [K]e,{R}e cho tung kieu PT co bien tudo khac nhau, PT khong co %
% bien tudo                                                               %  
% Khai bao dieukienbien: cac chuyen vi nut bang khong tuong ung voi cac   %      
% dieu kien: Mxy=Mx=My=Qx=Qy=0 trong {q}e cua {q}e=[A]e*{anfa}            %  
% Matran do cung nen tinh voi [B] cua PT khong co bien tudo               %  
% Thidu tinh tam vuong 6x6m                                               %  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%

clear all;
format long;     %dinh dang so
nel=144;         % So phan tu tam
nnel=4;          % So nut cho 1 phan tu
ndof=3;          % So bac tu do 1 nut
nnode=169;       %Tong so nut cua he
sdof=nnode*ndof; % Tong so bac tu do cua he
%-----------------------
% Khai bao so lieu tinh
%-----------------------
apt=0.5; % Chieu rong PT tam theo truc x [m]
bpt=0.5; % Chieu rong PT tam theo truc y [m]
%
% Dac trung tiet dien vat lieu tiet dien
el=2.65*10^6; %Mo dun dan hoi vat lieu beton [T/m2]
g=el/2;       % Mo dun truot vat lieu
nu=0.3;       % He so Poatxong
h=0.2;        % Chieu day ban [m]
hesocat=1;
hesonenk1=4*10^3; % heso nen vincle 
m=0.25;           % Khoi luong rieng vat lieu [T.(sec)2/m]
%
% Tai trong
q=-4; % Tai phan bo [T/m2]
%----------------
% Ket noi cac nut
%----------------
% nodes(phan tu, nut tam 1-4)= chi so nut HTD chung; 
for j=1:12
    for i=1:12
        nodes(i+(j-1)*12,1)=i+(j-1)*13;
        nodes(i+(j-1)*12,2)=i+1+(j-1)*13;
        nodes(i+(j-1)*12,3)=i+1+j*13;
        nodes(i+(j-1)*12,4)=i+j*13;
    end
end % ket thuc for j=1:12
%-----------------------------------------
% Toa do nut phan tu (He toa do dia phuong)
%-----------------------------------------
x(1)=0;   y(1)=0;
x(2)=apt; y(2)=0;
x(3)=apt; y(3)=bpt;
x(4)=0;   y(4)=bpt;
%------------------------------- 
%Dieu kien bien tu do theo chuvi
%-------------------------------
% Khaibaodieukienbien tudo: chuyenvinut tren bien bang khong voicackieu PT:
%   1)BienB1 y=0:            q1=Qy=0;  q2=My=0;  q4=Qy=0;   q5=My=0;
 
%   2)BienB2 y=bpt:          q7=Qy=0;  q8=My=0;  q10=Qy=0;  q11=My=0;
 
%   3)BienB3 x=0:            q1=Qx=0;  q3=Mx=0;; q10=Qx=0;  q12=Mx=0;
 
%   4)BienB4 x=apt:          q4=Qx=0;  q6=Mx=0;  q7=Qx=0;   q9=Mx=0;
 
%   5)BienB13 x=0 va y=0:    q1=Mxy=0; q2=My=0;  q3=Mx=0;   q4=Qy=0;
%                             q5=My=0;  q10=Qx=0; q12=Mx=0;
 
%   6)BienB23 x=0 va y=bpt:  q1=Qx=0;   q3=Mx=0;  q7=Qy=0;   q8=My=0;
%                             q10=Mxy=0; q11=My=0; q12=Mx=0;
 
%   7)BienB24 x=apt va y=0:  q4=Qx=0;  q6=Mx=0;  q7=Mxy=0;  q8=My=0;
%                             q9=Mx=0;  q10=Qy=0; q11=My=0;
 
%   8)BienB14 x=apt, y=apt:  q1=Qy=0;  q2=My=0;  q4=Mxy=0; q5=My=0;
%                             q6=Mx=0;  q7=Qx=0;  q9=Mx=0;
%   9) Bien Bo tuong ung voi phan tu khong co bien tu do
% Chi so qi (i=1-12) theo chi so chuyen vi nut trong HTD dia phuong
% So cvi nut bi rang buoc (bang khong)
i=0;
for j=1:13 % tu nut j=1-13
    if j==1 % Nut 1 
       i=i+1; 
       bcdof(i)=3*j-2;
       i=i+1; 
       bcdof(i)=3*j-1; 
       i=i+1; 
       bcdof(i)=3*j; 
    end % dktra 
    if j==13 % Nut 13
       i=i+1; 
       bcdof(i)=3*j-2;
       i=i+1; 
       bcdof(i)=3*j-1; 
       i=i+1; 
       bcdof(i)=3*j; 
    end 
    i=i+1; 
    bcdof(i)=3*j-2;
    i=i+1; 
    bcdof(i)=3*j-1; % Nut 2-12
end
for j=1:12   % Tu 14-157 
    if j==12 % nut 157
       i=i+1;  
       bcdof(i)=3*((j-1)*13+13)+1;
       i=i+1;  
       bcdof(i)=3*((j-1)*13+13)+2;
       i=i+1;  
       bcdof(i)=3*((j-1)*13+13)+3;
    end 
    i=i+1;  
    bcdof(i)=3*((j-1)*13+13)+1;
    i=i+1;  
    bcdof(i)=3*((j-1)*13+13)+3; % Nut 14-144
end
for j=1:12   % Tu 26-169 dung
    if j==12 %Nut 169
       i=i+1; 
       bcdof(i)=3*(j*13+13)-2;
       i=i+1; 
       bcdof(i)=3*(j*13+13)-1;
       i=i+1; 
       bcdof(i)=3*(j*13+13);
    end 
    i=i+1; 
    bcdof(i)=3*(j*13+13)-2;
    i=i+1; 
    bcdof(i)=3*(j*13+13);
end
for j=1:11 % Tu 158-168
    i=i+1; 
    bcdof(i)=3*(157+j)-2;
    i=i+1; 
    bcdof(i)=3*(157+j)-1;
end 
for inbc=1:i
    bcval(inbc)=0;
end
nbc=i;
% Ketthuc khai bao bien tudo theo chu vi
%
%-----------------------------------------------------------
% Kich thuoc va gan gia tri bang khong cho ma tran va vec to
%------------------------------------------------------------
indext=zeros(nnel*ndof,1);
ftq=zeros(12,1);     % vec to luc nut qui doi cua PT
ftqbo=zeros(12,1);   % vec to luc nut qui doi cua PT khong co bien tudo
ftqb1=zeros(12,1);   % vec to luc nut qui doi cua PT co bien tudo B1
ftqb2=zeros(12,1);   % vec to luc nut qui doi cua PT co bien tudo B2
ftqb3=zeros(12,1);   % vec to luc nut qui doi cua PT co bien tudo B3
ftqb4=zeros(12,1);   % vec to luc nut qui doi cua PT co bien tudo B4
ftqb13=zeros(12,1);  % vec to luc nut qui doi cua PT co bien tudo B13
ftqb23=zeros(12,1);  % vec to luc nut qui doi cua PT co bien tudo B23
ftqb24=zeros(12,1);  % vec to luc nut qui doi cua PT co bien tudo B24
ftqb14=zeros(12,1);  % vec to luc nut qui doi cua PT co bien tudo B14
%
ff=zeros(sdof,1);         % vec to luc nut cua he
matrana=zeros(12,12);     % [A] va [A]-1 theo (3.42)
matranabo=zeros(12,12);   % [A] va [A]-1 PT khong co bien tudo
matranab1=zeros(12,12);   % [A] va [A]-1 PT co bien tudo B1
matranab2=zeros(12,12);   % [A] va [A]-1 PT co bien tudo B2
matranab3=zeros(12,12);   % [A] va [A]-1 PT co bien tudo B3
matranab4=zeros(12,12);   % [A] va [A]-1 PT co bien tudo B4
matranab13=zeros(12,12);  % [A] va [A]-1 PT co bien tudo B13
matranab23=zeros(12,12);  % [A] va [A]-1 PT co bien tudo B23
matranab24=zeros(12,12);  % [A] va [A]-1 PT co bien tudo B24
matranab14=zeros(12,12);  % [A] va [A]-1 PT co bien tudo B14
matranc=zeros(3,12);      % Matran daohamrieng theo (3.47)
matranp=zeros(1,12);      % Matran [P(x,y)]
matranb=zeros(1,12);      % Matran [B] theo (3.44)
matranbbo=zeros(1,12);    % Matran [B] cua PT khong co bien tudo
matranbb1=zeros(1,12);    % Matran [B] cua PT co bien tudo B1
matranbb2=zeros(1,12);    % Matran [B] cua PT co bien tudo B2
matranbb3=zeros(1,12);    % Matran [B] cua PT co bien tudo B3
matranbb4=zeros(1,12);    % Matran [B] cua PT co bien tudo B4
matranbb13=zeros(1,12);   % Matran [B] cua PT co bien tudo B13
matranbb23=zeros(1,12);   % Matran [B] cua PT co bien tudo B23
matranbb24=zeros(1,12);   % Matran [B] cua PT co bien tudo B24
matranbb14=zeros(1,12);   % Matran [B] cua PT co bien tudo B14
kt=zeros(12,12);          % Matran docung PT tam
ktbo=zeros(12,12);        % Matran docung PT tam khong co bien tudo
ktb1=zeros(12,12);        % Matran docung PT tam co bien tudo B1
ktb2=zeros(12,12);        % Matran docung PT tam co bien tudo B2
ktb3=zeros(12,12);        % Matran docung PT tam co bien tudo B3
ktb4=zeros(12,12);        % Matran docung PT tam co bien tudo B4
ktb13=zeros(12,12);       % Matran docung PT tam co bien tudo B13
ktb23=zeros(12,12);       % Matran docung PT tam co bien tudo B23
ktb24=zeros(12,12);       % Matran docung PT tam co bien tudo B24
ktb14=zeros(12,12);       % Matran docung PT tam co bien tudo B14
ktnen1bo=zeros(12,12);    % Matran docung nen 01heso PT khongco bien tudo
kkt=zeros(sdof,sdof);     % Matran do cung he
% 
mmt=zeros(sdof,sdof);     % Matran khoiluong he
mt=zeros(12,12);          % Matran docung PT tam
mtbo=zeros(12,12);        % Matran docung PT tam khong co bien tudo
mtb1=zeros(12,12);        % Matran docung PT tam co bien tudo B1
mtb2=zeros(12,12);        % Matran docung PT tam co bien tudo B2
mtb3=zeros(12,12);        % Matran docung PT tam co bien tudo B3
mtb4=zeros(12,12);        % Matran docung PT tam co bien tudo B4
mtb13=zeros(12,12);       % Matran docung PT tam co bien tudo B13
mtb23=zeros(12,12);       % Matran docung PT tam co bien tudo B23
mtb24=zeros(12,12);       % Matran docung PT tam co bien tudo B24
mtb14=zeros(12,12);       % Matran docung PT tam co bien tudo B14
disp=zeros(sdof,1);       % vec to chvi nut he baitoan tinh
dispt=zeros(12,1);        % vec to cvi nut phan tu tam
%
noilucnut_tinh=zeros(3,4*nel); % Matran momen Mx,My,Mxy bai toan tinh
noilucnut_dong=zeros(3,4*nel); % Matran momen Mx,My,Mxy bai toan dong
noiluctg=zeros(3,1); % Matran trunggian momen Mx,My,Mxy bai toan tinh
%
%-------------------------------------------------------------------
% Xacdinh matranhangso [A]e cho cac kieu PT co va khong co bien tudo 
% ------------------------------------------------------------------
dp=el*h^3/12/(1-nu^2);% Do cung tru      
matranet=dp*[1  nu 0;...
             nu 1  0;...
             0  0  0.5*(1-nu)];% Ma tran [Et]
% Matran [A] Ptu khong co bien tudo - kyhieu [matranabo]
%-------------------------------------------------------
[matranabo] = bientudoB0(x, y);
% Matran [A] Ptu co bien tudo-B1 -kyhieu [matranab1]
%---------------------------------------------------
[matranab1] = bientudoB1(x,y,nu);
% Matran [A] Ptu co bien tudo-B2 -kyhieu [matranab2]
%---------------------------------------------------
[matranab2] = bientudoB2(x,y,nu);
% Matran [A] Ptu co bien tudo-B3 -kyhieu [matranab3]
%---------------------------------------------------
[matranab3] = bientudoB3(x,y,nu);
% Matran [A] Ptu co bien tudo-B4 -kyhieu [matranab4]
%---------------------------------------------------
[matranab4] = bientudoB4(x,y,nu);
% Matran [A] Ptu co bien tudo-B13 -kyhieu [matranab13]
%-----------------------------------------------------
[matranab13] = bientudoB13(x,y,nu);
% Matran [A] Ptu co bien tudo-B23 -kyhieu [matranab23]
%-----------------------------------------------------
[matranab23] = bientudoB23(x,y,nu);
% Matran [A] Ptu co bien tudo-B24 -kyhieu [matranab24]
%-----------------------------------------------------
[matranab24] = bientudoB24(x,y,nu);
% Matran [A] Ptu co bien tudo-B14 -kyhieu [matranab14]
%-----------------------------------------------------
[matranab14] = bientudoB14(x,y,nu);
%-----------------------------------------------------
% Tinh [K]e, {R}e - tinh tich phan bang cach doi bien 
% sang hetoado tunhien
%-----------------------------------------------------
% voi tam chu nhat doi bien sang HTD tu nhien: 
% x=0.5*apt*(r+1)   y=0.5*bpt*(s+1)  dxdy=0.25*apt*bpt*drds
for i=1:4 % tinh tich phan so 2x2 diem Gauss
    a=0.5773502691; 
    if i==1         
       r=-a;        
       s=-a;
    end
    if i==2
       r=-a;
       s=a;
    end
    if i==3
       r=a;
       s=-a;
    end
    if i==4
       r=a;
       s=a;
    end
    %[matranc]-matran daohamrieng bac 2 [Pxy] theo x2, y2 va xy(3.47)
    %[matranp]-matran [P(x,y)] theo (3.39)
    % 02 matran nay khong phu thuoc dieu kien bien
    matranc=[0 0 0 2 0 0 6*0.5*apt*(r+1) 2*0.5*bpt*(1+s) 0               0               6*0.5*apt*(r+1)*0.5*bpt*(1+s)       0;...
             0 0 0 0 0 2 0               0               2*0.5*apt*(r+1) 6*0.5*bpt*(1+s) 0                                   6*0.5*apt*(r+1)*bpt*(1+s);...
             0 0 0 0 2 0 0               4*0.5*apt*(r+1) 4*0.5*bpt*(1+s) 0               6*(0.5*apt*(r+1))^2                 6*(0.5*bpt*(1+s))^2];
    matranp=[1 0.5*apt*(r+1) 0.5*bpt*(1+s) (0.5*apt*(r+1))^2 0.5*apt*(r+1)*0.5*bpt*(1+s) (0.5*bpt*(1+s))^2 (0.5*apt*(r+1))^3 (0.5*apt*(r+1))^2*0.5*bpt*(1+s) 0.5*apt*(r+1)*(0.5*bpt*(1+s))^2 (0.5*bpt*(1+s))^3 (0.5*apt*(r+1))^3*0.5*bpt*(1+s) 0.5*apt*(r+1)*(0.5*bpt*(1+s))^3];
    %----------------------
    % PT khong co bien tudo
    matrandbo=matranc*matranabo; % Matran [D]e ngang 
    ktbo=ktbo+matrandbo'*matranet*matrandbo*0.25*apt*bpt; % [K]e 
    matranbbo=matranp*matranabo; % [B]e=[Pxy]*[A]-1  
    ftqbo=ftqbo+matranbbo'*q*0.25*apt*bpt; % {R}e
    ktnen1bo=ktnen1bo+hesonenk1*matranbbo'*matranbbo*0.25*apt*bpt;
    mtbo=mtbo+m*h*matranbbo'*matranbbo*0.25*apt*bpt; % [M]e
    %
    %----------------
    % PT bien tudo B1
    matrandb1=matranc*matranab1; % Matran [D]e ngang 
    ktb1=ktb1+matrandb1'*matranet*matrandb1*0.25*apt*bpt; % [K]e
    matranbb1=matranp*matranab1; % [B]e=[Pxy]*[A]-1
    ftqb1=ftqb1+matranbb1'*q*0.25*apt*bpt; % {R}e
    mtb1=mtb1+m*h*matranbb1'*matranbb1*0.25*apt*bpt; % [M]e
    %-----------------
    % PT bien tu do B2
    matrandb2=matranc*matranab2; % Matran [D]e ngang
    ktb2=ktb2+matrandb2'*matranet*matrandb2*0.25*apt*bpt; % [K]e
    matranbb2=matranp*matranab2; % [B]e=[Pxy]*[A]-1
    ftqb2=ftqb2+matranbb2'*q*0.25*apt*bpt; % {R}e
    mtb2=mtb2+m*h*matranbb2'*matranbb2*0.25*apt*bpt; % MT khoiluong PT B2
    %-----------------
    % PT bien tu do B3
    matrandb3=matranc*matranab3; % Matran [D]e ngang
    ktb3=ktb3+matrandb3'*matranet*matrandb3*0.25*apt*bpt; % [K]e
    matranbb3=matranp*matranab3; % [B]e=[Pxy]*[A]-1
    ftqb3=ftqb3+matranbb3'*q*0.25*apt*bpt; % {R}e
    mtb3=mtb3+m*h*matranbb3'*matranbb3*0.25*apt*bpt; % MT khoiluong PT B3
    %-----------------
    % PT bien tu do B4
    matrandb4=matranc*matranab4; % Matran [D]e ngang
    ktb4=ktb4+matrandb4'*matranet*matrandb4*0.25*apt*bpt; % [K]e
    matranbb4=matranp*matranab4; % [B]e=[Pxy]*[A]-1
    ftqb4=ftqb4+matranbb4'*q*0.25*apt*bpt; % {R}e
    mtb4=mtb4+m*h*matranbb4'*matranbb4*0.25*apt*bpt; % MT khoiluong PT B4
    %------------------
    % PT bien tu do B13
    matrandb13=matranc*matranab13; % Matran [D]e ngang
    ktb13=ktb13+matrandb13'*matranet*matrandb13*0.25*apt*bpt; % [K]e
    matranbb13=matranp*matranab13; % [B]e=[Pxy]*[A]-1
    ftqb13=ftqb13+matranbb13'*q*0.25*apt*bpt; % {R}e
    mtb13=mtb13+m*h*matranbb13'*matranbb13*0.25*apt*bpt; % MTkhoiluong PT B13
    %------------------
    % PT bien tu do B23
    matrandb23=matranc*matranab23; % Matran [D]e ngang
    ktb23=ktb23+matrandb23'*matranet*matrandb23*0.25*apt*bpt; % [K]e
    matranbb23=matranp*matranab23; % [B]e=[Pxy]*[A]-1
    ftqb23=ftqb23+matranbb23'*q*0.25*apt*bpt; % {R}e
    mtb23=mtb23+m*h*matranbb23'*matranbb23*0.25*apt*bpt; % MT khoiluong PT B23
    %------------------
    % PT bien tu do B24
    matrandb24=matranc*matranab24; % Matran [D]e ngang
    ktb24=ktb24+matrandb24'*matranet*matrandb24*0.25*apt*bpt; % [K]e
    matranbb24=matranp*matranab24; % [B]e=[Pxy]*[A]-1
    ftqb24=ftqb24+matranbb24'*q*0.25*apt*bpt; % {R}e
    mtb24=mtb24+m*h*matranbb24'*matranbb24*0.25*apt*bpt; % MT khoiluong PT B24
    %------------------
    % PT bien tu do B14
    matrandb14=matranc*matranab14; % Matran [D]e ngang
    ktb14=ktb14+matrandb14'*matranet*matrandb14*0.25*apt*bpt; % [K]e
    matranbb14=matranp*matranab14; % [B]e=[Pxy]*[A]-1
    ftqb14=ftqb14+matranbb14'*q*0.25*apt*bpt; % {R}e
    mtb14=mtb14+m*h*matranbb14'*matranbb14*0.25*apt*bpt; % MT khoiluong PT B14
end %for i=1:4 - Tinh cong don theo cac diem lay tich phan
% 
% Matran do cung nen (ktnen1bo) khong phan biet co bien tu do nhu doi voi PT tam
ktbo=ktbo+ktnen1bo;
ktb1=ktb1+ktnen1bo;
ktb2=ktb2+ktnen1bo;
ktb3=ktb3+ktnen1bo;
ktb4=ktb4+ktnen1bo;
ktb13=ktb13+ktnen1bo;
ktb23=ktb23+ktnen1bo;
ktb24=ktb24+ktnen1bo;
ktb14=ktb14+ktnen1bo;
%
% CHU Y: Matran docung nen nhu nhau voi cac kieu phantu, khong phu thuoc
%        dieu kien bien
%------------------------------------------------------------------------
% Ghep matran docung, matran khoiluong va vecto luc nut qui doi vao [kkt],[mmt], {ff} cua he
for iel=1:nel 
    % PT khong co bien tudo
    % Tinh [K]e=[Ktbo]e, {R}e={ftqbo}
    %--------------------------------
    kt=ktbo;
    ftq=ftqbo;
    mt=mtbo;
    % Bientudo B1 (PT 2 den 11)
    % Tinh [K]e=[Ktb1]e, {R}e={ftqb1}
    %--------------------------------
    if iel==2
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==3
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==4
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==5
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==6
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==7
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==8
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==9
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==10
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    if iel==11
       kt=ktb1;
       ftq=ftqb1;
       mt=mtb1;
    end
    % Bientudo B2 PT: 134 den 143
    % Tinh [K]e=[Ktb2]e, {R}e={ftqb2}
    %--------------------------------
    if iel==134
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==135
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==136
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==137
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==138
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==139
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==140
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==141
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==142
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    if iel==143
       kt=ktb2;
       ftq=ftqb2;
       mt=mtb2;
    end
    % Bientudo B3-PT:13,25,37,49,61,73,85,97,109,121
    % Tinh [K]e=[Ktb3]e, {R}e={ftqb3}
    %--------------------------------
    if iel==13
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==25
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==37
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==49
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==61
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==73
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==85
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==97
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==109
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    if iel==121
       kt=ktb3;
       ftq=ftqb3;
       mt=mtb3;
    end
    % Bientudo B4-PT:24,36,48,60,72,84,96,108,120,132
    % Tinh [K]e=[Ktb4]e, {R}e={ftqb4}
    %--------------------------------
    if iel==24
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==36
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==48
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==60
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==72
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==84
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==96
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==108
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==120
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    if iel==132
       kt=ktb4;
       ftq=ftqb4;
       mt=mtb4;
    end
    % Bientudo B13-PT:1
    % Tinh [K]e=[Ktb13]e, {R}e={ftqb13}
    %----------------------------------
    if iel==1
       kt=ktb13;
       ftq=ftqb13;
       mt=mtb13;
    end
    % Bientudo B23-PT133
    % Tinh [K]e=[Ktb23]e, {R}e={ftqb23}
    %----------------------------------
    if iel==133
       kt=ktb23;
       ftq=ftqb23;
       mt=mtb23;
    end
    % Bientudo B24-PT144
    % Tinh [K]e=[Ktb24]e, {R}e={ftqb24}
    %----------------------------------
    if iel==144
       kt=ktb24;
       ftq=ftqb24;
       mt=mtb24;
    end
    % Bientudo B14-PT12
    % Tinh [K]e=[Ktb14]e, {R}e={ftqb14}
    %----------------------------------
    if iel==12
       kt=ktb14;
       ftq=ftqb14;
       mt=mtb14;
    end
    for i=1:nnel    
        nd(i)=nodes(iel,i);
    end
    indext=tmhoacvinut(iel,nodes,nnel,ndof);% ma hoa cvi nut he
    mmt=ghepmtklt(mt,mmt,indext);
    kkt=ghepmtdct(kt,kkt,indext);
    edof=length(indext);
    for i=1:edof
        ii=indext(i);
        ff(ii)=ff(ii)+ftq(i);
    end
end % for iel=1:nel


%------------------------------------------------------------------------
% Giai baitoan tam chiu tai trong tinh
%------------------------------------------------------------------------
phantichtinh=1;
if phantichtinh==1 % Giai baitoan tinh hoc gan giatri: phantichtinh=1
   % Khu dieu kien bien
   [kkt,mmt,ff]=khudkbiendaodongt(kkt,mmt,ff,bcdof,bcval); 
   kkt=inv(kkt);
   disp=kkt*ff; 
   % --------------------------------------------------
   % Xac dinh noiluc phan tu tam Mx, My, Mxy cua PT tam
   for iel=1:nel
       %Khai bao loai phan tu
       matrana=matranabo; % [A]-1 MT nghich dao
       % Bientudo B1 (PT 2 den 11)
       %--------------------------
       if iel==2
          matrana=matranab1;
       end
       if iel==3
          matrana=matranab1;
       end
       if iel==4
          matrana=matranab1;
       end
       if iel==5
          matrana=matranab1;
       end
       if iel==6
          matrana=matranab1;
       end
       if iel==7
          matrana=matranab1;
       end
       if iel==8
          matrana=matranab1;
       end
       if iel==9
          matrana=matranab1;
       end
       if iel==10
          matrana=matranab1;
       end
       if iel==11
          matrana=matranab1;
       end
       % Bientudo B2 PT: 134 den 143
       %----------------------------
       if iel==134
          matrana=matranab2;
       end
       if iel==135
          matrana=matranab2;
       end
       if iel==136
          matrana=matranab2;
       end
       if iel==137
          matrana=matranab2;
       end
       if iel==138
          matrana=matranab2;
       end
       if iel==139
          matrana=matranab2;
       end
       if iel==140
          matrana=matranab2;
       end
       if iel==141
          matrana=matranab2;
       end
       if iel==142
          matrana=matranab2;
       end
       if iel==143
          matrana=matranab2;
       end
       % Bientudo B3-PT:13,25,37,49,61,73,85,97,109,121
       %-----------------------------------------------
       if iel==13
          matrana=matranab3;
       end
       if iel==25
          matrana=matranab3;
       end
       if iel==37
          matrana=matranab3;
       end
       if iel==49
          matrana=matranab3;
       end
       if iel==61
          matrana=matranab3;
       end
       if iel==73
          matrana=matranab3;
       end
       if iel==85
          matrana=matranab3;
       end
       if iel==97
          matrana=matranab3;
       end
       if iel==109
          matrana=matranab3;
       end
       if iel==121
          matrana=matranab3;
       end
       % Bientudo B4-PT:24,36,48,60,72,84,96,108,120,132
       %------------------------------------------------
       if iel==24
          matrana=matranab4;
       end
       if iel==36
          matrana=matranab4;
       end
       if iel==48
          matrana=matranab4;
       end
       if iel==60
          matrana=matranab4;
       end
       if iel==72
          matrana=matranab4;
       end
       if iel==84
          matrana=matranab4;
       end
       if iel==96
          matrana=matranab4;
       end
       if iel==108
          matrana=matranab4;
       end
       if iel==120
          matrana=matranab4;
       end
       if iel==132
          matrana=matranab4;
       end
       % Bientudo B13-PT:1
       %------------------
       if iel==1
          matrana=matranab13;
       end
       % Bientudo B23-PT133
       %-------------------
       if iel==133
          matrana=matranab23;
       end
       % Bientudo B24-PT144
       %--------------------
       if iel==144
          matrana=matranab24;
       end
       % Bientudo B14-PT12
       %------------------
       if iel==12
          matrana=matranab14;
       end
       %-----------------------------------------------------------------
       % Tinh chuyenvi w tai nut bien tudo tuong ung Qx hoac Qy hoac Mxy 
       % do khai bao dieukien bien bang khong. Khai bao:
       % iel=72; % Thidu tinh w tai nut 91
       % x=apt; % toa do (x,y) cua nut bien tinh w
       % y=bpt;
       % for i=1:4
       %    nd(i)=nodes(iel,i);
       %    dispt(3*i-2,1)=disp(3*nd(i)-2,1);
       %    dispt(3*i-1,1)=disp(3*nd(i)-1,1);
       %    dispt(3*i,1)=disp(3*nd(i),1);
       %end % for i=1:4   
       % khai bao matran [A]-1 tuong ung voi kieu phantu co nut bien
       %a=matranp*matranab4*dispt; % Giatri w tai nut bien tudo 91
       %-----------------------------------------------------------------
       % Tinh noiluc tai cac nut cua cac PT
       for iuon=1:4
           if iuon==1
              x=0;
              y=0;
           end
           if iuon==2
              x=apt;
              y=0;
           end
           if iuon==3
              x=apt;
              y=bpt;
           end
           if iuon==4
              x=0;
              y=bpt;
           end
           for i=1:4
               nd(i)=nodes(iel,i);
               dispt(3*i-2,1)=disp(3*nd(i)-2,1);
               dispt(3*i-1,1)=disp(3*nd(i)-1,1);
               dispt(3*i,1)=disp(3*nd(i),1);
           end % for i=1:4    
           matranpd=-[0 0 0 2 0 0 6*x 2*y 0 0 6*x*y 0;...
                      0 0 0 0 0 2 0 0 2*x 6*y 0 6*x*y;...
                      0 0 0 0 2 0 0 4*x 4*y 0 6*x^2 6*y^2];
           noiluctg=matranet*matranpd*matrana*dispt;
           noilucnut_tinh(:,4*(iel-1)+iuon)=noilucnut_tinh(:,4*(iel-1)+iuon)+noiluctg;      
       end %for iuon=1:4
   end % for iel=1:nel
end % ketthuc giai baitoan tinh if phantichtinh==1
%
%------------------------------------------------------------------------
% Giai baitoan dao dong - PP Newmark
%------------------------------------------------------------------------
pp_newmark=1;
if pp_newmark==1
   dieuhoa=0; %dieuhoa=1 Tinh voi tai trong dieu hoa gan giatri bang 1
   songnen=1; % songnen=1 Tinh voi tai trong songnen gan giatri bang 1
   % Khu dieu kien bien
   [kkt,mmt,ff]=khudkbiendaodongt(kkt,mmt,ff,bcdof,bcval);
   %----------------------------------------
   %
   % Tinh tan so daodong rieng va dang rieng
   fred=zeros(sdof-nbc,1);
   phi1=zeros(sdof-nbc,sdof-nbc);
   n=length(bcdof);
   phi=zeros(sdof,sdof);
   chiso=(1:sdof);
   chiso(bcdof)=[];
   kdr=kkt(chiso,chiso);
   mdr=mmt(chiso,chiso);
   [phi1,D]=eig(kdr,mdr);
   fred=sqrt(diag(D)); % tanso daodong rieng
   [fred,id]=sort(fred);% sapxep giatri tanso daodong rieng tu nho den lon
   %--------------------------------------------------- 
   %
   % Taitrong dieu hoa phanbo deu: q(x,y,t)=q*sin(omega*t)
   if dieuhoa==1
      ti=0;                % Thoi diem bat dau tinh
      tf=0.0314;           % thoi diem khao sat
      dt=0.0001;
      nt=fix((tf-ti)/dt);  % Tong so buoc thoi gian tinh
      omega=50;            % Tan so dao dong cuong buc 
      t=pi/100000;         % Thgian khao sat dang rieng
   end % if dieuhoa==1
   %
   %------------------------------------------------
   % Tai trong song nen phan bo deu cuongdo q theo qui luat:
   % 0 =< t =< teta1: f(t)=t/teta1 -  taitrong tang tuyen tinh
   % teta1=< t =< teta2: f(t)=1-t/teta2 - taitrong giam tuyen tinh
   % t > (teta1+teta2): f(t)=0 - taitrong bang khong
   if songnen==1
      teta1=0.1; % sec
      teta2=0.2; % sec
      teta3=0.01; % sec
      dt=0.0001; % buoc thoi gian (THAY DOI BUOC THOI GIAN TU 0.00001 GIAM 10 LAN XUONG CON 0.0001)
      
      % Khai bao tinh taitrong trong giai doan 
      giaidoan=1; % Tinh trong giaidoan taitrong tang gan gia tri 1
                  % Tinh trong giaidoan taitrong giam gan gia tri 2
                  % Tinh trong giaidoan taitrong bang khong gan gia tri 3
      %
      nt1=fix(teta1/dt);  % Tong so buoc thoi gian tinh giai doan 1
      nt2=fix(teta2/dt);  % Tong so buoc thoi gian tinh giai doan 2
      nt3=fix(teta3/dt);  % Tong so buoc thoi gian tinh giai doan 3
      if giaidoan==1             
         nt=nt1;
      end
      % tinh trong giai doan taitrong giam tuyen tinh
      if giaidoan==2 
         nt=nt1+nt2;
      end
         % tinh trong giai doan taitrong bang khong
      if giaidoan==3    
         nt=nt1+nt2+nt3;
      end   
   end % if songnen==1  
   disp_newmark=zeros(sdof,nt); % vectochuyenvinut cuahe- PP Newmark 
   acc=zeros(sdof,nt);          % acc la vec to gia toc cua he
   vel=zeros(sdof,nt);          % vel la vec to van toc cua he
   
   velo=zeros(sdof,1);  % vantoc ban dau
   acco=zeros(sdof,1);  % gia toc ban dau
   dispo=zeros(sdof,1); % chuyen vi ban dau
   khq=zeros(sdof,sdof);
   rhq=zeros(sdof,1);
   %  
   delta=0.5;
   afa=0.25;
   ao=1/afa/dt^2;
   a1=delta/afa/dt;
   a2=1/afa/dt;
   a3=(1/2/afa)-1;
   a4=(delta/afa)-1;
   a5=(dt/2)*((delta/afa)-2);
   a6=dt*(1-delta);
   a7=delta*dt;
   khq=kkt+ao*mmt; % matran docung hieu qua
   khq=inv(khq);
   %
   % Tinh vec to chuyenvi nut khi tam chiu tai trong dieu hoa
   if dieuhoa==1
      for it=1:nt
          if it==1
             rhq(:,1)=ff(:,1)*sin(omega*it*dt)+(mmt*(ao*dispo(:,1)+a2*velo(:,1)+a3*acco(:,1)));
             disp_newmark(:,it)=khq*rhq(:,1);
             acc(:,it)=ao*(disp_newmark(:,it)-dispo(:,1))-a2*velo(:,1)-a3*acco(:,1);
             vel(:,it)=velo(:,1)+a6*acco(:,1)+a7*acc(:,it);
          else
             rhq(:,1)=ff(:,1)*sin(omega*it*dt)+mmt*(ao*disp_newmark(:,it-1)+a2*vel(:,it-1)+a3*acc(:,it-1)); 
             disp_newmark(:,it)=khq*rhq(:,1);
             acc(:,it)=ao*(disp_newmark(:,it)-disp_newmark(:,it-1))-a2*vel(:,it-1)-a3*acc(:,it-1);
             vel(:,it)=vel(:,it-1)+a6*acc(:,it-1)+a7*acc(:,it);          
          end % Kthuc if it   
      end % for it=1:nt
   end %  if dieuhoa==1
   % 
      % Doan chuong trinh tinh w tai cac nut bien tudo chiu tai trong dieu hoa
      %------------------------------------------------------------------
      %iel=84; % phantu co nut bien tudo
      %x=apt; % toa do (x,y) cua nut bien tinh w
      %y=0;
      %for i=1:4
      %    nd(i)=nodes(iel,i);
      %    dispt(3*i-2,1)=disp_newmark(3*nd(i)-2,nt);
      %    dispt(3*i-1,1)=disp_newmark(3*nd(i)-1,nt);
      %    dispt(3*i,1)=disp_newmark(3*nd(i),nt);
      %end % for i=1:4  
       % khai bao matran [A]-1 tuong ung voi kieu phantu co nut bien
      %chuyenvi_nuttudo=matranp*matranab4*dispt; % Giatri w tai nut bien tudo
  %---------------------------------------------------------
  % Tinh vec to chuyenvi nut khi tam chiu tai trong song nen    
 if songnen==1
    % Tinh voi taitrong tang tuyen tinh - giai doan 1
    %------------------------------------------------
    if giaidoan==1
       for it=1:nt1 % giai doan tai trong tang tuyen tinh
           if it==1
              rhq(:,1)=ff(:,1)*(it*dt/teta1)+(mmt*(ao*dispo(:,1)+a2*velo(:,1)+a3*acco(:,1)));
              disp_newmark(:,it)=khq*rhq(:,1);
              acc(:,it)=ao*(disp_newmark(:,it)-dispo(:,1))-a2*velo(:,1)-a3*acco(:,1);
              vel(:,it)=velo(:,1)+a6*acco(:,1)+a7*acc(:,it);
           else
             rhq(:,1)=ff(:,1)*(it*dt/teta1)+mmt*(ao*disp_newmark(:,it-1)+a2*vel(:,it-1)+a3*acc(:,it-1)); 
             disp_newmark(:,it)=khq*rhq(:,1);
             acc(:,it)=ao*(disp_newmark(:,it)-disp_newmark(:,it-1))-a2*vel(:,it-1)-a3*acc(:,it-1);
             vel(:,it)=vel(:,it-1)+a6*acc(:,it-1)+a7*acc(:,it);   
           end % Kthuc if it==1   
       end % for it=1:nt1
    end % if giaidoan==1
    %
    % Tinh voi taitrong giam tuyen tinh - giai doan 2
    %------------------------------------------------
    if giaidoan==2
       for it=1:nt1 % giai doan tai trong tang tuyen tinh
           if it==1
rhq(:,1)=ff(:,1)*(it*dt/teta1)+(mmt*(ao*dispo(:,1)+a2*velo(:,1)+a3*acco(:,1)));
              disp_newmark(:,it)=khq*rhq(:,1);
              acc(:,it)=ao*(disp_newmark(:,it)-dispo(:,1))-a2*velo(:,1)-a3*acco(:,1);
              vel(:,it)=velo(:,1)+a6*acco(:,1)+a7*acc(:,it);
           else
             rhq(:,1)=ff(:,1)*(it*dt/teta1)+mmt*(ao*disp_newmark(:,it-1)+a2*vel(:,it-1)+a3*acc(:,it-1)); 
             disp_newmark(:,it)=khq*rhq(:,1);
             acc(:,it)=ao*(disp_newmark(:,it)-disp_newmark(:,it-1))-a2*vel(:,it-1)-a3*acc(:,it-1);
             vel(:,it)=vel(:,it-1)+a6*acc(:,it-1)+a7*acc(:,it);   
           end % Kthuc if it==1   
       end % for it=1:nt
       %
       % Thong so ban dau tinh cho giai doan 2
       dispo(:,1)=disp_newmark(:,nt1);
       velo(:,1)=vel(:,nt1);
       acco(:,1)=acc(:,nt1);
       %
       for it=1:nt2 % giai doan tai trong giam tuyen tinh
           if it==1
              rhq(:,1)=ff(:,1)*(1-it*dt/teta2)+(mmt*(ao*dispo(:,1)+a2*velo(:,1)+a3*acco(:,1)));
              disp_newmark(:,nt1+it)=khq*rhq(:,1);
              acc(:,nt1+it)=ao*(disp_newmark(:,nt1+it)-dispo(:,1))-a2*velo(:,1)-a3*acco(:,1);
              vel(:,nt1+it)=velo(:,1)+a6*acco(:,1)+a7*acc(:,nt1+it);
           else
              rhq(:,1)=ff(:,1)*(1-it*dt/teta2)+mmt*(ao*disp_newmark(:,(nt1+it)-1)+a2*vel(:,(nt1+it)-1)+a3*acc(:,(nt1+it)-1));
              disp_newmark(:,nt1+it)=khq*rhq(:,1);
              acc(:,nt1+it)=ao*(disp_newmark(:,nt1+it)-disp_newmark(:,(nt1+it)-1))-a2*vel(:,(nt1+it)-1)-a3*acc(:,(nt1+it)-1);
              vel(:,nt1+it)=vel(:,(nt1+it)-1)+a6*acc(:,(nt1+it)-1)+a7*acc(:,(nt1+it));          
           end % Kthuc if it==1   
       end % for it=1:nt2
    end % if giaidoan==2
    %
    % Tinh trong giaidoan 3: taitrong bang khong
    if giaidoan==3
       for it=1:nt1 % Tinh giai doan tai trong tang tuyen tinh
           if it==1
              rhq(:,1)=ff(:,1)*(it*dt/teta1)+(mmt*(ao*dispo(:,1)+a2*velo(:,1)+a3*acco(:,1)));
              disp_newmark(:,it)=khq*rhq(:,1);
              acc(:,it)=ao*(disp_newmark(:,it)-dispo(:,1))-a2*velo(:,1)-a3*acco(:,1);
              vel(:,it)=velo(:,1)+a6*acco(:,1)+a7*acc(:,it);
           else  
              rhq(:,1)=ff(:,1)*(it*dt/teta1)+mmt*(ao*disp_newmark(:,it-1)+a2*vel(:,it-1)+a3*acc(:,it-1));
              disp_newmark(:,it)=khq*rhq(:,1);
              acc(:,it)=ao*(disp_newmark(:,it)-disp_newmark(:,it-1))-a2*vel(:,it-1)-a3*acc(:,it-1);
              vel(:,it)=vel(:,it-1)+a6*acc(:,it-1)+a7*acc(:,it);          
           end % Kthuc if it==1   
       end % for it=1:nt
       % 
       % Thong so ban dau tinh cho giai doan 2
       dispo(:,1)=disp_newmark(:,nt1);
       velo(:,1)=vel(:,nt1);
       acco(:,1)=acc(:,nt1);
       %
       for it=1:nt2 % giai doan tai trong giam tuyen tinh
           if it==1
              rhq(:,1)=ff(:,1)*(1-it*dt/teta2)+(mmt*(ao*dispo(:,1)+a2*velo(:,1)+a3*acco(:,1)));
              disp_newmark(:,nt1+it)=khq*rhq(:,1);
              acc(:,nt1+it)=ao*(disp_newmark(:,nt1+it)-dispo(:,1))-a2*velo(:,1)-a3*acco(:,1);
              vel(:,nt1+it)=velo(:,1)+a6*acco(:,1)+a7*acc(:,nt1+it);
           else
              rhq(:,1)=ff(:,1)*(1-it*dt/teta2)+mmt*(ao*disp_newmark(:,(nt1+it)-1)+a2*vel(:,(nt1+it)-1)+a3*acc(:,(nt1+it)-1));
              disp_newmark(:,nt1+it)=khq*rhq(:,1);
              acc(:,nt1+it)=ao*(disp_newmark(:,nt1+it)-disp_newmark(:,(nt1+it)-1))-a2*vel(:,(nt1+it)-1)-a3*acc(:,(nt1+it)-1);
              vel(:,nt1+it)=vel(:,(nt1+it)-1)+a6*acc(:,(nt1+it)-1)+a7*acc(:,(nt1+it));          
           end % Kthuc if it==1   
       end % for it=1:nt2
       %
       % Thong so ban dau tinh cho giai doan 3
       dispo(:,1)=disp_newmark(:,nt1+nt2);
       velo(:,1)=vel(:,nt1+nt2);
       acco(:,1)=acc(:,nt1+nt2);
       %
       for it=1:nt3 % giai doan tai trong bang khong
           if it==1
              rhq(:,1)=(mmt*(ao*dispo(:,1)+a2*velo(:,1)+a3*acco(:,1)));
              disp_newmark(:,nt1+nt2+it)=khq*rhq(:,1);
              acc(:,nt1+nt2+it)=ao*(disp_newmark(:,nt1+nt2+it)-dispo(:,1))-a2*velo(:,1)-a3*acco(:,1);
              vel(:,nt1+nt2+it)=velo(:,1)+a6*acco(:,1)+a7*acc(:,nt1+nt2+it);
           else
              rhq(:,1)=mmt*(ao*disp_newmark(:,(nt1+nt2+it)-1)+a2*vel(:,(nt1+nt2+it)-1)+a3*acc(:,(nt1+nt2+it)-1));
              disp_newmark(:,nt1+nt2+it)=khq*rhq(:,1);
              acc(:,nt1+nt2+it)=ao*(disp_newmark(:,nt1+nt2+it)-disp_newmark(:,(nt1+nt2+it)-1))-a2*vel(:,(nt1+nt2+it)-1)-a3*acc(:,(nt1+nt2+it)-1);
              vel(:,nt1+nt2+it)=vel(:,(nt1+nt2+it)-1)+a6*acc(:,(nt1+nt2+it)-1)+a7*acc(:,(nt1+nt2+it));          
           end % Kthuc if it==1   
       end % for it=1:nt2
    end % if giaidoan==3 
    % Doan chuong trinh tinh w tai cac nut bien tudo chiu song nen- thi du tai nut 91
    %--------------------------------------------------------------------
    iel=84; % phantu co nut bien tudo
    x=apt; % toa do (x,y) cua nut bien tinh w
    y=0;
    for i=1:4
        nd(i)=nodes(iel,i);
        dispt(3*i-2,1)=disp_newmark(3*nd(i)-2,nt);
        dispt(3*i-1,1)=disp_newmark(3*nd(i)-1,nt);
        dispt(3*i,1)=disp_newmark(3*nd(i),nt);
    end % for i=1:4  
       % khai bao matran [A]-1 tuong ung voi kieu phantu co nut bien
    chuyenvi_nuttudo1=matranp*matranab4*dispt; % Giatri w tai nut bien tudo
    %
    iel=72; % phantu co nut bien tudo
    x=apt; % toa do (x,y) cua nut bien tinh w
    y=bpt;
    for i=1:4
        nd(i)=nodes(iel,i);
        dispt(3*i-2,1)=disp_newmark(3*nd(i)-2,nt);
        dispt(3*i-1,1)=disp_newmark(3*nd(i)-1,nt);
        dispt(3*i,1)=disp_newmark(3*nd(i),nt);
    end % for i=1:4  
       % khai bao matran [A]-1 tuong ung voi kieu phantu co nut bien
    chuyenvi_nuttudo2=matranp*matranab4*dispt; % Giatri w tai nut bien tudo
 end % if songnen==1
 % 
 % Tinh noi luc
 %-------------
 for iel=1:nel
     %Khai bao loai phan tu
     matrana=matranabo; % [A]-1 MT nghich dao
     % Bientudo B1 (PT 2 den 11)
     %--------------------------------
     if iel==2
        matrana=matranab1;
       end
       if iel==3
          matrana=matranab1;
       end
       if iel==4
          matrana=matranab1;
       end
       if iel==5
          matrana=matranab1;
       end
       if iel==6
          matrana=matranab1;
       end
       if iel==7
          matrana=matranab1;
       end
       if iel==8
          matrana=matranab1;
       end
       if iel==9
          matrana=matranab1;
       end
       if iel==10
          matrana=matranab1;
       end
       if iel==11
          matrana=matranab1;
       end
       % Bientudo B2 PT: 134 den 143
       %----------------------------
       if iel==134
          matrana=matranab2;
       end
       if iel==135
          matrana=matranab2;
       end
       if iel==136
          matrana=matranab2;
       end
       if iel==137
          matrana=matranab2;
       end
       if iel==138
          matrana=matranab2;
       end
       if iel==139
          matrana=matranab2;
       end
       if iel==140
          matrana=matranab2;
       end
       if iel==141
          matrana=matranab2;
       end
       if iel==142
          matrana=matranab2;
       end
       if iel==143
          matrana=matranab2;
       end
       % Bientudo B3-PT:13,25,37,49,61,73,85,97,109,121
       %-----------------------------------------------
       if iel==13
          matrana=matranab3;
       end
       if iel==25
          matrana=matranab3;
       end
       if iel==37
          matrana=matranab3;
       end
       if iel==49
          matrana=matranab3;
       end
       if iel==61
          matrana=matranab3;
       end
       if iel==73
          matrana=matranab3;
       end
       if iel==85
          matrana=matranab3;
       end
       if iel==97
          matrana=matranab3;
       end
       if iel==109
          matrana=matranab3;
       end
       if iel==121
          matrana=matranab3;
       end
       % Bientudo B4-PT:24,36,48,60,72,84,96,108,120,132
       %------------------------------------------------
       if iel==24
          matrana=matranab4;
       end
       if iel==36
          matrana=matranab4;
       end
       if iel==48
          matrana=matranab4;
       end
       if iel==60
          matrana=matranab4;
       end
       if iel==72
          matrana=matranab4;
       end
       if iel==84
          matrana=matranab4;
       end
       if iel==96
          matrana=matranab4;
       end
       if iel==108
          matrana=matranab4;
       end
       if iel==120
          matrana=matranab4;
       end
       if iel==132
          matrana=matranab4;
       end
       % Bientudo B13-PT:1
       %------------------
       if iel==1
          matrana=matranab13;
       end
       % Bientudo B23-PT133
       %-------------------
       if iel==133
          matrana=matranab23;
       end
       % Bientudo B24-PT144
       %--------------------
       if iel==144
          matrana=matranab24;
       end
       % Bientudo B14-PT12
       %------------------
       if iel==12
          matrana=matranab14;
       end
       % Tinh noiluc tai cac nut cua cac PT
       for iuon=1:4
           if iuon==1
              x=0;
              y=0;
           end
           if iuon==2
              x=apt;
              y=0;
           end
           if iuon==3
              x=apt;
              y=bpt;
           end
           if iuon==4
              x=0;
              y=bpt;
           end
           % CHU y: nt tinh voi taitrong dieu hoa
           % khi tinh voi song nen can khai bao gia tri ntinh
           if songnen==1
              ntinh=1000; % tinh voi tai song nen THAY DOI TU 15000 XUONG 1000
                        % ntinh <= nt1 - giaidoan 1
                        % nt1 <=ntinh<= nt1+nt2 - giaidoan 2
                        % nt1+nt2 <=ntinh<= nt1+nt2+nt3 - giaidoan 3
              nt=ntinh;
         end
           
           for i=1:4
               nd(i)=nodes(iel,i);
               dispt(3*i-2,1)=disp_newmark(3*nd(i)-2,nt);
               dispt(3*i-1,1)=disp_newmark(3*nd(i)-1,nt);
               dispt(3*i,1)=disp_newmark(3*nd(i),nt);
           end % for i=1:4 
           matranpd=-[0 0 0 2 0 0 6*x 2*y 0 0 6*x*y 0;...
                      0 0 0 0 0 2 0 0 2*x 6*y 0 6*x*y;...
                      0 0 0 0 2 0 0 4*x 4*y 0 6*x^2 6*y^2];
           noiluctg=matranet*matranpd*matrana*dispt;
           noilucnut_dong(:,4*(iel-1)+iuon)=noilucnut_dong(:,4*(iel-1)+iuon)+noiluctg;      
       end %for iuon=1:4
   end % for iel=1:nel
end % Kthuc phantich dong - PP Newmark
%-------------------------------------
