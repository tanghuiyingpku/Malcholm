function  isGrow = FracturePropagation(CritTheta,CritK,KI,KII)
global TipStates TipStatesInv AllEle_global ConnList_global ;
global KIChf KICnf  PresF_global  nActTip 
global ArrivalNF InitialAperture Mat;
global RecordK ;
global Kindex;
global NFgrow epsKI;
NFgrow =0;
KIChf = 1;
KICnf = 1;
KIChf_tip = 1;
% Relative Value
HFinit = 0;%InitialAperture;
isGrow = 0;
%Cross with Extisting Nf- path judgement
TipStates0 = TipStates;
nTip0 = nActTip;
for itip = 1 : nTip0
    %Open HF at the end of NF
    nf = TipStatesInv(itip);
    if AllEle_global(nf,10) > 2.2  && TipStates0(nf) < 998 &&  ConnList_global(nf,3)*ConnList_global(nf,4) < -1e-5
        iParent = ConnList_global(nf,1);
       % if ArrivalNF(iParent) > 0.1 && TipStates0(nf) < 990
        if  TipStates0(nf) < 990
            if  (abs(CritK(TipStates0(nf))- KIChf_tip) < epsKI || CritK(TipStates0(nf))- KIChf_tip >= epsKI) %&& KII(TipStates0(nf)) > 5E-1
           %     DrawDisplacement_new();
                RecordK(Kindex) = CritK(TipStates0(nf))/KIChf_tip;
                Kindex =Kindex + 1;
                fprintf('Open %dth Tip\n',TipStates0(nf));
                isGrow = 1;
                disp('Left Tip Grow...');
                Pinit = PresF_global(nf)/5;%Pres(nf) - abs(Q(1))*AllEle(1,7)*12/TipOpening(nf)^3/FracHeight;
                Dinit = -HFinit; %-e_global(nf);%-HFinit;
                DynamicGrowth_nfTip(0,Pinit,Dinit,nf,CritTheta(TipStates0(nf)));
                %Dtemp = CalcDD(nFracture,Fractures,nAllEle,AllEle,ConnList,Mat,Ptemp,EleType);
                continue;
                
            end
        end
    end
    %kink grow
    
    
    kinkG = 0;
    if ((AllEle_global(nf,10) < 1.1 || abs(AllEle_global(nf,10)-6) < 1e-6 ) && TipStates0(nf) < 998)
        if  abs(CritK(TipStates0(nf))- KIChf)/KIChf < epsKI|| CritK(TipStates0(nf))- KIChf >= 0.1
            ncc = ConnList_global(nf,2);
            for ibc = 1 : ncc
                nb = ConnList_global(nf,2+ibc);
                if nb > 0.1
                    if AllEle_global(nb,10) > 6
                        kinkG = 1;
                    end
                end
            end
            if kinkG > 0.1
                fprintf('Open %dth Tip\n',TipStates0(nf));
                RecordK(Kindex) = CritK(TipStates0(nf))/KIChf;
                Kindex =Kindex + 1;
                isGrow = 1;
                disp('Left Tip Grow...');
                Pinit = PresF_global(nf);%
                Dinit =  -HFinit;%DD_global(MaxEle+nf)/2;%
                DynamicGrowth_kink(Pinit,Dinit,nf);
                continue;
            end
        end
    end
        
    %Open HF
    %modify
    if ((AllEle_global(nf,10) < 1.1 || abs(AllEle_global(nf,10)-6) < 1e-6 ) && TipStates0(nf) < 998) 
        if  abs(CritK(TipStates0(nf))- KIChf)/KIChf < epsKI|| CritK(TipStates0(nf))- KIChf >= 0.1
            fprintf('Open %dth Tip\n',TipStates0(nf));
            RecordK(Kindex) = CritK(TipStates0(nf))/KIChf;
            Kindex =Kindex + 1;
            isGrow = 1;
            disp('Left Tip Grow...');
            Pinit = PresF_global(nf);%
            Dinit =  -HFinit;%DD_global(MaxEle+nf)/2;%
            DynamicGrowth(0,Pinit,Dinit,nf,CritTheta(TipStates0(nf)));
            continue;
        end
    end
    
    OpenNF = 0;
    Deviate_Coeff = -0.1;
    Along_Coeff =0;
    %Open HF at the middle of NF
    if (AllEle_global(nf,10) > 2.1  && TipStates0(nf) < 998) &&  ConnList_global(nf,3)*ConnList_global(nf,4)>1e-5
        if  abs(CritK(TipStates0(nf))- KIChf)/KIChf < epsKI || CritK(TipStates0(nf))- KIChf >= epsKI && KI(TipStates0(nf)) > 1E-6
            % if Connect with another HF not propagate
            Deviate_Coeff = 0;%CritK(TipStates0(nf))/KIChf;
            OpenNF = 0;
        end
    end
    % Open NF
    if (AllEle_global(nf,10) > 2.1  && TipStates0(nf) < 998) &&  ConnList_global(nf,3)*ConnList_global(nf,4)>1e-5
        iParent = ConnList_global(nf,1);
        % Along NF
        if ArrivalNF(iParent) > 0.1
            Gnf = (KI(TipStates0(nf))/Mat.K1NF)^2 +  (KII(TipStates0(nf))/Mat.K2NF)^2;
            if  abs(Gnf- KICnf) < epsKI ||Gnf -  KICnf>= epsKI
                Along_Coeff = Gnf/ KICnf;
                OpenNF = 1;
            end
%             if  abs(abs(KI2(TipStates0(nf)))- KICnf) < epsKI ||KI2(TipStates0(nf))- KICnf>= epsKI
%                 Along_Coeff = abs(KI2(TipStates0(nf)))/ KICnf;
%                 OpenNF = 1;
%             end
        end
    end
    if OpenNF>0.1
        NFgrow = 1;
        if Along_Coeff > Deviate_Coeff
            RecordK(Kindex) = Along_Coeff;
            Kindex =Kindex + 1;
            fprintf('Open %dth Tip\n',TipStates0(nf));
            isGrow = 1;
            disp('Left Tip Grow...');
            Pinit = PresF_global(nf)*0.5;
            Dinit = 0;
            DynamicGrowthNF(Pinit,Dinit,nf);
            continue;
        else
            %if Along_Coeff > 1.001
            continue;
            tipcor = Tipcoordinate(nf,:);
            if abs(AllEle_global(nf,1) - tipcor(1)) < 1e-6 && abs(AllEle_global(nf,2) - tipcor(2)) < 1e-6
                tiphf = [AllEle_global(nf,3) AllEle_global(nf,4)];
            else
                tiphf = [AllEle_global(nf,1) AllEle_global(nf,2)];
            end
            if Deviate_Coeff < 1e-10
                continue;
            end
            RecordK(Kindex) = Deviate_Coeff;
            Kindex =Kindex + 1;
            fprintf('Open %dth Tip Step Over\n',TipStates0(nf));
            isGrow = 1;
            disp('Left Tip Grow...');
            Pinit = PresF_global(nf);%
            Dinit = 0;%
            DynamicGrowth_Stepover(Pinit,Dinit,nf,CritTheta(TipStates0(nf)));
            continue;
            % end
        end
    end
end
isGrow2 = isGrow;
% The boundary NF elements are actived
isGrow = NaturalEleActivation2(isGrow2);
isAct = 0;%NaturalEleActivation();
if isAct > 0.1
    isGrow = 1;
end
end
