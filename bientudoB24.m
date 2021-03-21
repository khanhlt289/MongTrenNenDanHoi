function [matranab24] = bientudoB24(x,y,nu);
% Matran [A] Ptu co bien tudo-B24 -kyhieu [matranab24]
% Ham chuyen vi xap xi bang da thuc Pascal
% Ngay kiem tra: 3-7-2002
% Nguoi kiem tra: Luu Truong Khanh
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
matranab24=[1 0  0 0    0     0      0         0         0         0           0              0;...
            0 0  1 0    0     0      0         0         0         0           0              0;...
            0 -1 0 0    0     0      0         0         0         0           0              0;...
            0 0  0 0    0     0      6         0         2*(2-nu)  0           0              0;...
            0 0  1 0    x(2)  0      0         x(2)^2    0         0           x(2)^3         0;...
            0 0  0 2    0     2*nu   6*x(2)    0         2*nu*x(2) 0           0              0;...
            0 0  0 0    1     0      0         2*x(3)    2*y(3)    0           3*x(3)^2       3*y(3)^2;...
            0 0  0 2*nu 0     2      6*nu*x(3) 2*nu*y(3) 2*x(3)    6*y(3)      6*nu*x(3)*y(3) 6*x(3)*y(3);...
            0 0  0 2    0     2*nu   6*x(3)    2*y(3)    2*nu*x(3) 6*nu*y(3)   6*x(3)*y(3)    6*nu*x(3)*y(3);...
            0 0  0 0    0     0      0         2*(2-nu)  0         6           0              0;...
            0 0  0 2*nu 0     2      0         2*nu*y(4) 0         6*y(4)      0              0;...
            0 -1 0 0    -y(4) 0      0         0         -y(4)^2   0           0              -y(4)^3];
matranab24=inv(matranab24);