function MassBalanceCheck(nAllEle,AllEle,DD,Pres,t,nwell,well,Vinit,Mat)
Vinj = 0;
for i = 1 : nwell
    Vinj = Vinj+well{i}.q*t;
end
Vinj = Vinj + Vinit;
disp('Injection Fluid Volume is:');
fprintf('    %f   \n',Vinj);
Vstore = CalcStoreV(nAllEle,AllEle,Pres*1e6,DD,Mat);
disp('Stored Fluid Volume is:');
fprintf('    %f   \n',Vstore);
fprintf(' Relative Mass Error is %f\n',abs(Vstore-Vinj)/Vinj);
fprintf(' Absolute Mass Error is %f\n',Vstore-Vinj);
end