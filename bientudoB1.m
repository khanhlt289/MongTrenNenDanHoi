function [matranab1] = bientudoB1(x,y,nu);
% Matran [A] Ptu co bien tudo-B1 -kyhieu [matranab1]
% Ham chuyen vi xap xi bang da thuc Pascal
% Ngay kiem tra: 3-7-2002
% Nguoi kiem tra: Luu Truong Khanh
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
matranab1=[0 0    0    0       0         0      0           2*(2-nu)     0             6        0              0;...
           0 0    0    2*nu    0         2      0           0            0             0        0              0;...
           0 -1   0    0       0         0      0           0            0             0        0              0;...
           0 0    0    0       0         0      0           2*(2-nu)     0             6        6*(2-nu)*x(2)  6*x(2);...
           0 0    0    2*nu    0         2      6*nu*x(2)   0            2*x(2)        0        0              0;...
           0 -1   0    -2*x(2) 0         0      -3*x(2)^2   0            0             0        0              0;...
           1 x(3) y(3) x(3)^2  x(3)*y(3) y(3)^2 x(3)^3      x(3)^2*y(3)  x(3)*y(3)^2   y(3)^3   x(3)^3*y(3)    x(3)*y(3)^3;...
           0 0    1    0       x(3)      2*y(3) 0           x(3)^2       2*x(3)*y(3)   3*y(3)^2 x(3)^3         3*x(3)*y(3)^2;...
           0 -1   0    -2*x(3) -y(3)     0      -3*x(3)^2   -2*x(3)*y(3) -y(3)^2       0        -3*x(3)^2*y(3) -y(3)^3;...
           1 0    y(4) 0       0         y(4)^2 0           0            0             y(4)^3   0              0;...
           0 0    1    0       0         2*y(4) 0           0            0             3*y(4)^2 0              0;...
           0 -1   0    0       -y(4)     0      0           0            -y(4)^2       0        0              -y(4)^3];      
matranab1=inv(matranab1);
%
