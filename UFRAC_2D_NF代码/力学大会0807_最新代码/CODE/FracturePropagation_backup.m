function  isGrow = FracturePropagation_backup(CritTheta,CritK,KI1,isGrow0)
global TipStates TipStatesInv AllEle_global ConnList_global ;
global nAct MaxEle nTip
global KIChf KICnf Mat PresF_global  DD_global nActTip Tipcoordinate
global ArrivalNF InitialAperture;
global RecordK;
global Kindex;
global NFgrow;
NFgrow =0;
KIChf = 1;
KIChf_tip = 1;
% Relative Value
HFinit = InitialAperture;
isGrow = 0;
%Cross with Extisting Nf- path judgement
TipStates0 = TipStates;
nAct0 = nAct;
isOpenNf = 0;
nTip0 = nActTip;
for itip = 1 : nTip0
    %Open HF at the end of NF
    nf = TipStatesInv(itip);
    if AllEle_global(nf,10) > 2.2  && TipStates0(nf) < 998 &&  ConnList_global(nf,3)*ConnList_global(nf,4) < -1e-5 
        iParent = ConnList_global(nf,1);
        if ArrivalNF(iParent) > 0.1 && TipStates0(nf) < 990
            if  abs(CritK(TipStates0(nf))- KIChf_tip) < 0.1 || CritK(TipStates0(nf))- KIChf_tip >= 0.1
                RecordK(Kindex) = CritK(TipStates0(nf))/KIChf_tip;
                if RecordK(Kindex)>5
                    f = 1;
                end
                Kindex =Kindex + 1;
                d = AllEle_global(nf,7)/2;
                K1 = -Mat.G*1e3/4/(1-Mat.miu)*(2*pi/d)^0.5*DD_global(MaxEle+nf);
                if K1 > 1e-3
                    fprintf('Open %dth Tip\n',TipStates0(nf));
                    isGrow = 1;
                    isBranch = 1;
                    disp('Left Tip Grow...');
                    Pinit = PresF_global(nf)/5;%Pres(nf) - abs(Q(1))*AllEle(1,7)*12/TipOpening(nf)^3/FracHeight;
                    Dinit = -HFinit; %-e_global(nf);%-HFinit;
                    DynamicGrowth_nfTip(0,Pinit,Dinit,nf,CritTheta(TipStates0(nf)));
                    %Dtemp = CalcDD(nFracture,Fractures,nAllEle,AllEle,ConnList,Mat,Ptemp,EleType);
                    continue;
                end
                
            end
        end
    end
    %Open HF
    if (AllEle_global(nf,10) < 1.1  && TipStates0(nf) < 998)
        if  abs(CritK(TipStates0(nf))- KIChf)/KIChf <0.1 || CritK(TipStates0(nf))- KIChf >= 0.1
            fprintf('Open %dth Tip\n',TipStates0(nf));
            RecordK(Kindex) = CritK(TipStates0(nf))/KIChf;
            Kindex =Kindex + 1;
            isGrow = 1;
            disp('Left Tip Grow...');
            Pinit = PresF_global(nf)*0.5;%
            Dinit =  -HFinit;%DD_global(MaxEle+nf)/2;%
            DynamicGrowth(0,Pinit,Dinit,nf,CritTheta(TipStates0(nf)));
            %Dtemp = CalcDD(nFracture,Fractures,nAllEle,AllEle,ConnList,Mat,Ptemp,EleType);
            continue;
        end
    end
    OpenNF = 0;
    Deviate_Coeff = -0.1;
    Along_Coeff =0;
    %Open HF at the middle of NF
    if (AllEle_global(nf,10) > 2.1  && TipStates0(nf) < 998) &&  ConnList_global(nf,3)*ConnList_global(nf,4)>1e-5
        if  abs(CritK(TipStates0(nf))- KIChf)/KIChf <0.1 || CritK(TipStates0(nf))- KIChf >= 0.1
            % if Connect with another HF not propagate
            Deviate_Coeff = CritK(TipStates0(nf))/KIChf;
            OpenNF = 1;
        end
    end  
    % Open NF
    if (AllEle_global(nf,10) > 2.1  && TipStates0(nf) < 998) &&  ConnList_global(nf,3)*ConnList_global(nf,4)>1e-5
        iParent = ConnList_global(nf,1);
        if ArrivalNF(iParent) > 0.1
            if  abs(KI1(TipStates0(nf))- KICnf) < 0.1 ||KI1(TipStates0(nf))- KICnf>= 0.1
                Along_Coeff = KI1(TipStates0(nf))/ KICnf;
                OpenNF = 1;
            end
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
            Pinit = PresF_global(nf);
            Dinit = 0;
            DynamicGrowthNF(Pinit,Dinit,nf);
            continue;
        else
            %continue;
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
        end
    end 
end
% The boundary NF elements are actived?
isAct = 0;%NaturalEleActivation();
if isAct > 0.1
    isGrow = 1;
end
end
