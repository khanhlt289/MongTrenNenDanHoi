function [matranab3] = bientudoB3(x,y,nu);
% Ham chuyen vi xap xi bang da thuc Pascal
% Ngay kiem tra: 3-7-2002
% Nguoi kiem tra: Luu Truong Khanh
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
matranab3=[0 0    0    0       0         0      6         0            2*(2-nu)    0         0              0;...
           0 0    1    0       0         0      0         0            0           0         0              0;...
           0 0    0    2       0         2*nu   0         0            0           0         0              0;...
           1 x(2) 0    x(2)^2  0         0      x(2)^3    0            0           0         0              0;...
           0 0    1    0       x(2)      0      0         x(2)^2       0           0         x(2)^3         0;...
           0 -1   0    -2*x(2) 0         0      -3*x(2)^2 0            0           0         0              0;...
           1 x(3) y(3) x(3)^2  x(3)*y(3) y(3)^2 x(3)^3    x(3)^2*y(3)  x(3)*y(3)^2 y(3)^3    x(3)^3*y(3)    x(3)*y(3)^3;...
           0 0    1    0       x(3)      2*y(3) 0         x(3)^2       2*x(3)*y(3) 3*y(3)^2  x(3)^3         3*x(3)*y(3)^2;...
           0 -1   0    -2*x(3) -y(3)     0      -3*x(3)^2 -2*x(3)*y(3) -y(3)^2     0         -3*x(3)^2*y(3) -y(3)^3;...
           0 0    0    0       0         0      6         0            2*(2-nu)    0         6*y(4)         6*(2-nu)*y(4);...
           0 0    1    0       0         2*y(4) 0         0            0           3*y(4)^2  0              0;...
           0 0    0    2       0         2*nu   0         2*y(4)       0           6*nu*y(4) 0              0]; 
matranab3=inv(matranab3);
%