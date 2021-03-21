%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
% Khu dieu kien bien bai toan dao dong la rut gon ma tran do cung tam,    %  
% ma tran khoi luong va vec to tai quy nut co nhieu phan tu khac 0 hon    %  
% so voi o ma tran dia phuong                                             %
% Nguoi thuc hien: Luu Truong Khanh                                       %
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
function [kkt,mmt,ff]=khudkbiendaodongt(kkt,mmt,ff,bcdof,bcval)           %  
n=length(bcdof);                                                          %      
sdof=size(kkt);                                                           %  
for i=1:n;                                                                %  
    c=bcdof(i);                                                           %  
    for j=1:sdof;                                                         %  
        kkt(c,j)=0;                                                       %  
        kkt(j,c)=0;                                                       %  
        mmt(c,j)=0;                                                       %  
        mmt(j,c)=0;                                                       %  
    end;                                                                  %  
    kkt(c,c)=1;                                                           %  
    mmt(c,c)=1;                                                           %  
    ff(c)=bcval(i);                                                       %  
end;                                                                      %  