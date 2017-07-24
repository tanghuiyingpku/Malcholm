function InitialGlobalVaria(paraFile)
global EleType;
global Mat;
global Fluid;
global CM_global
global GM_global
global CM_marker
global GrowthN;
global MaxEle nTStep;
global UseTipE;
global alpha;
global beta;
global UseHeight;
global InitialAperture; %Inital Aperture of natural fractures
global TipStates TipStatesInv Tipcoordinate;
global stepL DT_init TEND;
global KIChf KICnf
global isMechActive_global
global PicScale isRefineTip NumofRefine;
global CpF_global XfF_global Hslurry_global Hbanking_global Tau_global;
global PresF_global DD_global e_global AllEle_global ConnList_global sigmaN_global;
global G_inj;
global epsKI;
epsKI = 0.1;
G_inj = 0;
fid = fopen(paraFile,'r');
temp = fscanf(fid,'%s',1);
while ~feof(fid)
    if strcmp(temp,'CONTROL') == 1
        while strcmp(temp,'END') < 0.1
            temp = fscanf(fid,'%s',1);
            if strcmp(temp,'ELETYPE') == 1
                EleType = fscanf(fid,'%d',1);
            end
            if strcmp(temp,'USEHEIGHT') == 1
                UseHeight = fscanf(fid,'%d',1);
            end
            if strcmp(temp,'USETIP') == 1
                UseTipE = fscanf(fid,'%d',1);
            end
            if strcmp(temp,'MAXElE') == 1
                MaxEle = fscanf(fid,'%d',1);
            end
            if strcmp(temp,'MAXTIMESTEP') == 1
                nTStep = fscanf(fid,'%d',1);
            end
            if strcmp(temp,'DT') == 1
                DT_init = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'TIME') == 1
                TEND = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'PICSCALE') == 1
                PicScale = fscanf(fid,'%f %f %f %f',4);
            end
            if strcmp(temp,'REFINETIP') == 1
                isRefineTip = fscanf(fid,'%d',1);
                NumofRefine = fscanf(fid,'%d',1);
            end
        end
    end
    if strcmp(temp,'RESERVOIR') == 1
        while strcmp(temp,'END') < 0.1
             temp = fscanf(fid,'%s',1);
            if strcmp(temp,'SXX') == 1
                Mat.Sxx = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'SYY') == 1
                Mat.Syy = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'SXY') == 1
                Mat.Sxy = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'Sv') == 1
                Mat.Sv = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'SMINDIR') == 1
                Mat.SminDir = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'PP') == 1
                Mat.Pp = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'E') == 1
                Mat.E = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'PR') == 1
                Mat.miu = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'HEIGHT') == 1
                Mat.h = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'TENS') == 1
                Mat.tens = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'PORO') == 1
                Mat.poro = fscanf(fid,'%f',1);
            end
        end
      %  Mat.E = Mat.E/(1-Mat.miu^2);
        Mat.G = Mat.E/2/(1+Mat.miu);
    end
    if strcmp(temp,'STIMULATION') == 1
        while strcmp(temp,'END') < 0.1
             temp = fscanf(fid,'%s',1);
            if strcmp(temp,'PROPPANT') == 1
                Fluid.nprop = fscanf(fid,'%d',1);
                CpF_global = zeros(MaxEle,Fluid.nprop);
                for ii = 1 : Fluid.nprop
                    Fluid.prop{ii}.name = fscanf(fid,'%s',1);
                    Fluid.prop{ii}.diam = fscanf(fid,'%f',1);
                    Fluid.prop{ii}.dens = fscanf(fid,'%f',1);
                end
            end
            if strcmp(temp,'FLUID') == 1
                Fluid.nfluid = fscanf(fid,'%d',1);
                XfF_global = zeros(MaxEle,Fluid.nfluid);
                for ii = 1 : Fluid.nfluid
                    Fluid.fluid{ii}.name = fscanf(fid,'%s',1);
                    Fluid.fluid{ii}.dens = fscanf(fid,'%f',1);
                    Fluid.fluid{ii}.compr = fscanf(fid,'%f',1);
                    Fluid.fluid{ii}.viscotype = fscanf(fid,'%s',1);
                    Fluid.fluid{ii}.refp = fscanf(fid,'%f',1);
                    Fluid.fluid{ii}.refp = Fluid.fluid{ii}.refp/1e6;
                end
            end
            if strcmp(temp,'RHEOLOGY') == 1
                for ii = 1 : Fluid.nfluid
                    name = fscanf(fid,'%s',1);
                    for jj = 1 : Fluid.nfluid
                        if strcmp(name, Fluid.fluid{jj}.name) == 1
                            Fluid.fluid{jj}.visco = fscanf(fid,'%f',1);
                            break;
                        end
                    end
                end
            end
            if strcmp(temp,'LEAKOFF') == 1
                Fluid.spurtslop = fscanf(fid,'%f',1);
                Fluid.spurtdeltap = fscanf(fid,'%f',1);
                Fluid.spurtv = fscanf(fid,'%f',1);
                for ii = 1 : Fluid.nfluid
                    name = fscanf(fid,'%s',1);
                    for jj = 1 : Fluid.nfluid
                        if strcmp(name, Fluid.fluid{jj}.name) == 1
                            Fluid.fluid{jj}.isleakoff = fscanf(fid,'%f',1);
                            break;
                        end
                    end
                end
            end
        end
    end
    if strcmp(temp,'FRAC') == 1
        while strcmp(temp,'END') < 0.1
             temp = fscanf(fid,'%s',1);
            if strcmp(temp,'FRIC') == 1
                Mat.fric = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'FULLPEN') == 1
                Mat.fullpen = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'KICHF') == 1
                KIChf = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'KICNF') == 1
                KICnf = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'STEPL') == 1
                stepL = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'INITIALWIDTH_HF') == 1
                InitialAperture = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'GROWTHN') == 1
                GrowthN = fscanf(fid,'%f',1);
            end
            alpha = 1;
            beta = 2.3;
        end
    end
    if strcmp(temp,'ANISOTROPY') == 1
        while strcmp(temp,'END') < 0.1
            temp = fscanf(fid,'%s',1);
            if strcmp(temp,'K1HF') == 1
                Mat.K1HF = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'K2HF') == 1
                Mat.K2HF = fscanf(fid,'%f',1);
            end
             if strcmp(temp,'K1NF') == 1
                Mat.K1NF = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'K2NF') == 1
                Mat.K2NF = fscanf(fid,'%f',1);
            end
            if strcmp(temp,'BETA') == 1
                Mat.Kbeta = fscanf(fid,'%f',1);
            end
        end
      %  fprintf('**** Calculate Rotation Angle Matrix ****');
      %  FindMaxTheta2(Mat.Kmax,Mat.Kmin)
    end
    temp = fscanf(fid,'%s',1);
end
% Open Space 

Tipcoordinate = zeros(MaxEle,2) - 999;
TipStatesInv = zeros(MaxEle,1);
isMechActive_global = zeros(MaxEle,1)+1;
CM_global = zeros(MaxEle*2,MaxEle*2)-1;
GM_global = zeros(MaxEle,MaxEle);
Hslurry_global = zeros(MaxEle,Fluid.nprop)+Mat.h;
Hbanking_global =zeros(MaxEle,Fluid.nprop);
TipStates = zeros(MaxEle,1);
Tau_global = zeros(MaxEle,1);
sigmaN_global  = zeros(MaxEle,1);
% 几何位置的标记；TIP 2 新增单元 3 与TIP相邻 1 与TIP的邻居相邻 0  内部 -1
CM_marker = zeros(MaxEle,1)-1;
PresF_global = zeros(MaxEle,1);
DD_global = zeros(MaxEle*2,1);
e_global  = zeros(MaxEle,1);
AllEle_global = zeros(MaxEle,11);
ConnList_global = zeros(MaxEle,8);

