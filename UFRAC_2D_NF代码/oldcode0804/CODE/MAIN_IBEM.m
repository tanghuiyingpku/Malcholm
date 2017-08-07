function   MAIN_IBEM()
tic;
ClearData();
Casename = 'test';
clc;
close('all');
clear global ;
disp('Delete Old Pictures');
delete([pwd,'\*.png']);
fclose('all');
global FILEPATH  FILEPATH_main;
%Casename = 'multiperf';
MAINPATH = '..\';
FILEPATH = [MAINPATH,'DATA\',Casename,'\Fig\'];
ResultPath = [MAINPATH,'DATA\',Casename,'\Data\'];
mkdir(ResultPath);
delete([ResultPath,'*.in']);
mkdir(FILEPATH);
FILEPATH_main =  [MAINPATH,'DATA\',Casename,'\'];
fracFile =[FILEPATH_main, 'Input_frac.in'];
wellInput = [FILEPATH_main,'Input_well.in'];
paraFile = [FILEPATH_main,'Input_para.in'];
delete([FILEPATH,'\*.fig']);
delete([FILEPATH,'\*.bmp']);

disp('%*******************************************************');
disp('% *Indirect Boundary Element Method by Huiying Tang @2014');
disp('% *加入高度修正G');
disp('% *加入Joints连接');
disp('% *Joints包括三类情况：2 闭合胶结 3 开启 4 滑移');
disp('% *OBJ:验证curve fracture');
disp('%*******************************************************');
%% ************** Input Control Parameters *******************
disp('  ');
disp('Define Global Variables');
disp('  ');
%% *************** global variables ***********************
global nTStep ;
global HasInter;
global  nAct
global well nwell;
global DT_init TEND;
global ActWellGrid;
global RecordK;
global Kindex NFgrow;
global nActOld MassBL;
global GrowNumber einit_global;
global nAllEle_global   AllEle_global   ConnList_global

GrowNumber = 4;
nActOld = 0;
Kindex =1;
RecordK = zeros(1000,1);
MassBL = zeros(1000,1);
einit_global = zeros(9000,1);

%% Define global variable
InitialGlobalVaria(paraFile);
%% Read well file and write fracture files
ReadWellFile(wellInput);
disp('  ');
fprintf('There are %d wells in total\n',nwell);
disp('  ');
% Whether NF enters HF
%% Read properties files 
disp('Reading Input File');
ReadInputFile(fracFile);
disp('  ');
%% Initialization
disp('Find Corresponding well grid index');
FindWellIndex();
%Initialize Fracture Part - active numbers
Initialization_couple();

dt0 = DT_init;
timelist = zeros(nTStep+1,1);
timelist(1) = dt0;

nt = 1;
% IterFile = 'Iteration.dat';
% delete(IterFile);
isGrow = 0;
CurT = 0;
nG = well{1}.nGrid;
ActWellGrid = zeros(nG*nwell,1);
dt = dt0;
nGrow = 0;
NFgrow= 0;

IsInit = 0;
GrowStep = 0;
tic;
tt = 0;
hasR = 0;
dt0 = dt;
HasInter=0;
DrawDs(0,0);

for i = 1 : nAllEle_global
    if ConnList_global(i,2) == 3
        hold on;
        plot(AllEle_global(i,8), AllEle_global(i,9),'r*');
    end
    if ConnList_global(i,2) == 4
        hold on;
        plot(AllEle_global(i,8), AllEle_global(i,9),'g*');
    end
end
%% Start Simulation
while CurT < TEND
    %Well Part
%     if mod(nt,5) ==0
%         SearchRadius();
%     end
    for ii =1 : nwell
        % Solve Each well seperately
        IsInit = CheckFracInit(1e8,ii,CurT);
    end
    if nAct < 0.1
        nActOld = 0;
        continue;
    else
        %% Solve  Equations
        if IsInit > 0.1
            %Update DD according to new Pressure
            %CalcDD();
            isGrow = 1;
        end
        if nGrow < 0.1
            dt = dt0/10;
        end
        if hasR > 0.1
            tt = tt+ 1;
        end
        if GrowStep < 0
            dt =  FluidSolidFullCouple_Scaling(CurT,dt,isGrow);
        else
            dt =  FluidSolidFullCouple_Scaling_well(CurT,dt,isGrow,nt);
        end
        
        % Calculate stress intensity factor
        [KI1,KI2] = StressIntensF();
        [CritTheta,CritG] = FindCritDir_Anisotropy(KI1,KI2);
        CurT = CurT+dt;
        if dt > dt0*10
            dt = dt0*10;
        end
        if dt < dt0/10
            dt = dt0/10;
        end
        isGrow = 0;
        nActOld =  nAct;
        if HasInter > 0.1
            dt = dt0;
            disp('Judge Intersection Behavior with Energy Release Rate');
            [hasRV,hasR,nf]  = HasReachCritical(CritG,KI1);
            if hasR > 0.1
                GrowthPathJudge_GIGII2(nf,CritTheta,KI1,KI2);
                close all;
                isGrow = 1;
                dt = dt / 5;
            end
        end
        GrowStep = GrowStep+1;
        %% Fracture propagation
        if isGrow < 0.1
            isGrow = FracturePropagation(CritTheta,CritG,KI1,KI2);
        end
        if isGrow > 0.1
            title(['nt = ',num2str(nt)]);
            UpdateTau(nActOld, nAct,CurT);
            nGrow = nGrow + 1;
        end
    end
    %% Save well data
    saveWellData(CurT);
    err = MassBalance_global(CurT);
    MassBL(nt) = err;
    timelist(nt) = CurT;
    nt = nt + 1;
    isEnd();
    if CurT > 201
        f = 1;
    end
end
%% Final timestep
toc;
saveData();
figure
plot(RecordK(1:Kindex),'.');
fclose all;
% DrawDn_t(0,0)

end
