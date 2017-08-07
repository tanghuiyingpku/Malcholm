function  UpdateMass()
% The boundary condition for fluid flow in no flow bounday
% In this function Shear Displacements has nothing to do with fluid flow so
% is not taken into direct calculation
% The elements with no fluid flow is not taken into calculation
%% Get global variables
global  IndexInv nAct   
global   Mat  ;
global  MaxEle ;
global DD_global PresF_global AllEle_global ConnList_global;
global   e_global  einit_global ;

clear isMechActive0;

%% Open Space for variables
P_pre = zeros(nAct,1);
AllEle = zeros(nAct,11);
ConnList = zeros(nAct,8);
DD_pre = zeros(nAct*2,1);
epsd = 1e-8;


%Scaling characteristic value
Maxtau = 12*1e6;
einit = zeros(nAct,1);
%% Assign values
for i = 1 : nAct
    DD_pre(i) = DD_global(IndexInv(i));
    DD_pre(i+nAct) = DD_global(MaxEle+IndexInv(i));
    P_pre(i) = PresF_global(IndexInv(i))*1e6;
    AllEle(i,:) = AllEle_global(IndexInv(i),:);
    ConnList(i,:) = ConnList_global(IndexInv(i),:);
    %alpha = theta + 90;
    einit(i) = e_global(IndexInv(i));
end

%friction coefficient
fric = tand(Mat.fric);
[~,~,~,e0,~] = CheckStates(einit,nAct, AllEle, epsd,DD_pre,P_pre,Mat, fric, Maxtau);
% isMechActive0 = isMechActive0*0+1;
e = e0;
for i =  1 : nAct
    if abs(einit(i) - e(i)) > 1e-8 %%&& AllEle(i,10) > 1.1 && einit(i) < 1e-8
        einit_global(IndexInv(i)) = e(i);
    end
end

%% give value back to global variabls
for i = 1 : nAct
    e_global(IndexInv(i)) =e(i);
end

clear e0 e P_pre  DD_pre DD DDcur;
clear  InsituS   AllEle   ;

end


