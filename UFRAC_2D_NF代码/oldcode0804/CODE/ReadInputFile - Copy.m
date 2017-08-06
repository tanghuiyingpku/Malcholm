%% 读取输入文件，同时进行裂缝网格剖分和赋值(加密）
function  ReadInputFile(inputFile)
global TipStates;
global CM_marker;
global ArrivalNF;
global nTip Tipcoordinate;
global Fractures;
global nFracture nHF_global;
global AllEle_global ConnList_global nAllEle_global ;
nTip = 0;

minEleL = 3*1e-4;
ConnList_global(:,2) = 2;
fid = fopen(inputFile,'r');
% 裂缝个数
fid = ReadComment(fid);
nFracture = fscanf(fid,'%d',1);
fprintf('There are %d Fractures\n',nFracture);
%建立Fracture的结构体
%包括{Type, 属性值（4*1）
Fracture = cell(nFracture,1);
FracLoc = zeros(nFracture,4);
% 裂缝类型
fid = ReadComment(fid);
fprintf('The Fracture Types are :\n');
ArrivalNF = zeros(nFracture,1)*(-1);
nHF_global = 0;
for i = 1 : nFracture
    Fracture{i}.Type =  fscanf(fid,'%d',1);
    if Fracture{i}.Type == 1
        disp('Hydraulic Fracture');
        nHF_global = nHF_global + 1;
    else
        ArrivalNF(i) = 0;
        disp('Natural Fracture/Joints');
    end
end
% 每个裂缝的初始平均单元长度
fid = ReadComment(fid);
for i = 1 : nFracture
    Fracture{i}.EleLAve = fscanf(fid,'%f',1);
end
% 裂缝位置
fid = ReadComment(fid);
for i = 1 : nFracture
    FracLoc(i,1) = fscanf(fid,'%f',1);
    FracLoc(i,2) = fscanf(fid,'%f',1);
    FracLoc(i,3) = fscanf(fid,'%f',1);
    FracLoc(i,4) = fscanf(fid,'%f',1);
    if FracLoc(i,3) < FracLoc(i,1)
        temp = FracLoc(i,3);
        FracLoc(i,3)=FracLoc(i,1);
        FracLoc(i,1) = temp;
        temp = FracLoc(i,4);
        FracLoc(i,4)=FracLoc(i,2);
        FracLoc(i,2) = temp;
    end
    vec = [(FracLoc(i,4)-FracLoc(i,2)),(FracLoc(i,3)-FracLoc(i,1))];
    % 约定角度在第一和第四象限
    Fracture{i}.length = CalculateDis(FracLoc(i,1:2),FracLoc(i,3:4));
    if vec(1)*vec(2) >=  0
        Fracture{i}.sint = abs(vec(1))/sqrt(vec(1)^2+ vec(2)^2);
        Fracture{i}.cost = abs(vec(2))/sqrt(vec(1)^2+ vec(2)^2);
    else
        Fracture{i}.sint = -abs(vec(1))/sqrt(vec(1)^2+ vec(2)^2);
        Fracture{i}.cost = abs(vec(2))/sqrt(vec(1)^2+ vec(2)^2);
    end
    Fracture{i}.Loc = FracLoc(i,:);
end
% for i = 1 : nFracture
%     plot([FracLoc(i,1) FracLoc(i,3)],[FracLoc(i,2) FracLoc(i,4)],'m');
%     hold on;
% end
%裂缝属性
fid = ReadComment(fid);
for i = 1 : nFracture
    % 压裂裂缝
    if Fracture{i}.Type == 1
        Fracture{i}.Ss = fscanf(fid,'%f',1);
        Fracture{i}.Sn = fscanf(fid,'%f',1);
    end
    if Fracture{i}.Type == 7
        Fracture{i}.Ss = fscanf(fid,'%f',1);
        Fracture{i}.Sn = fscanf(fid,'%f',1);
    end
    % 天然裂缝
    if Fracture{i}.Type ==2
        Fracture{i}.Ks = fscanf(fid,'%f',1);
        Fracture{i}.Kn = fscanf(fid,'%f',1);
        Fracture{i}.Fang = fscanf(fid,'%f',1);
        Fracture{i}.Coh = fscanf(fid,'%f',1);
    end
end
%读取完毕
fclose(fid);
nAllEle_global = 0;
%寻找裂缝间的相交关系并按顺序记录相交的网格编号
countEle = 1;
countTip = 1;
for i = 1 : nFracture
    startn = countEle;
    CM_marker(startn) = 2;
    CM_marker(startn+1) = 1;
    CM_marker(startn+2) = 0;
    %Preset the space
    Point = zeros(50,2);
    NeighbNum = zeros(50,1);
    nCrossPoint = 0;
    L0 = [FracLoc(i,1),FracLoc(i,2),FracLoc(i,3),FracLoc(i,4)];
    for j = 1 : nFracture
        if j ~= i
            Li = [FracLoc(j,1),FracLoc(j,2),FracLoc(j,3),FracLoc(j,4)];
            % isCross: whether cross CrossPoint where to cross
            [isCross,CrossPoint] = FindCrossPoint(L0,Li);
            % Adjust length
            % If 2 Fracs are too close, bound them together, change one
            % end to be the cross point
            lengthI1 = sqrt((Li(3)-Li(1))^2 + (Li(4)-Li(2))^2);
            lengthI2 = sqrt((L0(3)-L0(1))^2 + (L0(4)-L0(2))^2);
            lengthI = min(lengthI1,lengthI2);
            
            epss = 0.002;
            if abs(CrossPoint(1) - L0(1))^2 + abs(CrossPoint(2) - L0(2))^2 < 1*epss
                FracLoc(i,1) = CrossPoint(1);
                FracLoc(i,2) = CrossPoint(2);
                
              %  isCross=1;
            end
            if abs(CrossPoint(1) - L0(3))^2 + abs(CrossPoint(2) - L0(4))^2 < 1*epss
                FracLoc(i,3) = CrossPoint(1);
                FracLoc(i,4) = CrossPoint(2);
               
              %   isCross=1;
            end
            if abs(CrossPoint(1) - Li(1))^2 + abs(CrossPoint(2) - Li(2))^2 < 1*epss
                FracLoc(j,1)= CrossPoint(1);
                FracLoc(j,2) = CrossPoint(2);
               
                % isCross=1;
            end
            if abs(CrossPoint(1) - Li(3))^2 + abs(CrossPoint(2) - Li(4))^2 < 1*epss
                FracLoc(j,3) = CrossPoint(1);
                FracLoc(j,4) = CrossPoint(2);
               
              %  isCross=1;
            end
            [isCross,CrossPoint] = FindCrossPoint(FracLoc(i,1:4),FracLoc(j,1:4));
            if isCross > 0.1
                nCrossPoint = nCrossPoint+1;
                Point(nCrossPoint,:) = CrossPoint;
                NeighbNum(nCrossPoint) = j;
            end
        end
    end
    Point  = Point(1:nCrossPoint,:);
    NeighbNum = NeighbNum(1:nCrossPoint);
    [~,sequen] = sort(Point(:,1));
    Point = Point(sequen,:);
    NeighbNum = NeighbNum(sequen);
    Fracture{i}.nConnectF = nCrossPoint;
    Fracture{i}.NeighbF = NeighbNum;
    Fracture{i}.ConnectP = Point;
    %计算每个小段的单元个数
    aveLinit = Fracture{i}.EleLAve;
    nInit = round(Fracture{i}.length/aveLinit);
    AveL = zeros(nCrossPoint+1,1)+aveLinit;
    nEle = zeros(nCrossPoint+1,1);
    Li = zeros(nCrossPoint+1,1);
    %标记每个交点左侧的element
    index = zeros(nCrossPoint,1);
    %给单元编号
    Fracture{i}.Tip(1) = nAllEle_global+1;
    if Fracture{i}.Type < 1.1
        TipStates(nAllEle_global+1) = countTip;
        nTip = nTip + 1;
        countTip = countTip+1;
    end
    for ii = 1 : nCrossPoint
        if ii < 1.1
            ConnList_global(nAllEle_global+1,3) = -1;
            ConnList_global(nAllEle_global+1,4) = nAllEle_global+2;
            Li(ii) = CalculateDis(Point(ii,:),FracLoc(i,1:2));
            Tipcoordinate(nAllEle_global+1,1:2) = FracLoc(i,1:2);
        else
            Li(ii) = CalculateDis(Point(ii,:),Point(ii-1,:));
            if Li(ii) > 1e-6
                ConnList_global(nAllEle_global+1,3) = nAllEle_global;
                ConnList_global(nAllEle_global+1,4) = nAllEle_global+2;
            end
        end
        % AveL(ii) = Li(ii)/Fracture{i}.length*aveLinit;
        if AveL(ii) < 1e-6
            %两个交点重合
            nEle(ii)=0;
        else
            if AveL(ii) < minEleL
                AveL(ii) = minEleL;
            end
            nEle(ii) = ceil(Li(ii)/AveL(ii));
        end
        if nEle(ii) > 0.1
            %parent fracture
            ConnList_global(nAllEle_global+1:nAllEle_global+nEle(ii),1) = i;
        end
        for jj = 2 : nEle(ii)
            ConnList_global(nAllEle_global+jj,3) = nAllEle_global+jj-1;
            ConnList_global(nAllEle_global+jj,4) = nAllEle_global+jj+1;
        end
        nAllEle_global = nAllEle_global + nEle(ii);
        index(ii) = nAllEle_global;
    end
    if nAllEle_global > 420
        f = 1;
    end
    if nCrossPoint > 0.1
        Li(nCrossPoint+1) = CalculateDis(Point(nCrossPoint,:),FracLoc(i,3:4));
        % AveL(nCrossPoint+1) = Li(nCrossPoint+1)/Fracture{i}.length*aveLinit;
        if AveL(nCrossPoint+1) < minEleL
            AveL(nCrossPoint+1) = minEleL;
        end
        nEle(nCrossPoint+1) = ceil(Li(nCrossPoint+1)/AveL(nCrossPoint+1));
        for jj = 1 : nEle(nCrossPoint+1)-1
            ConnList_global(nAllEle_global+jj,3) = nAllEle_global+jj-1;
            ConnList_global(nAllEle_global+jj,4) = nAllEle_global+jj+1;
        end
        ConnList_global(nAllEle_global+1:nAllEle_global+nEle(nCrossPoint+1),1) = i;
        jj = jj + 1;
        if nEle(nCrossPoint+1) < 1.1
            jj = 1;
        end
        ConnList_global(nAllEle_global+jj,3) = nAllEle_global+jj-1;
        ConnList_global(nAllEle_global+jj,4) = -1;
        Tipcoordinate(nAllEle_global+1,1:2) = FracLoc(i,3:4);
        nAllEle_global = nAllEle_global + nEle(nCrossPoint+1);
        Fracture{i}.Li = Li;
        Fracture{i}.eleNum = nEle;
        Fracture{i}.index = index;
    else
        %no crosspoint???????
        Fracture{i}.index = nAllEle_global+1;
        nEle = nInit;
        ConnList_global(nAllEle_global+1:nAllEle_global+nEle,1) = i;
        ConnList_global(nAllEle_global+1,3) = -1;
        Tipcoordinate(nAllEle_global+1,1:2) = FracLoc(i,1:2);
        ConnList_global(nAllEle_global+1,4) = nAllEle_global+2;
        for jj = 2 : nEle-1
            ConnList_global(nAllEle_global+jj,3) = nAllEle_global+jj-1;
            ConnList_global(nAllEle_global+jj,4) = nAllEle_global+jj+1;
        end
        Fracture{i}.eleNum = nEle;
        nAllEle_global = nAllEle_global + nEle;
        ConnList_global(nAllEle_global,3) = nAllEle_global-1;
        ConnList_global(nAllEle_global,4) = -1;
        Tipcoordinate(nAllEle_global,1:2) = FracLoc(i,3:4);
    end
    Fracture{i}.Tip(2) = nAllEle_global;
    if Fracture{i}.Type < 1.1
        TipStates(nAllEle_global) = countTip;
        nTip = nTip + 1;
        countTip = countTip+1;
    end
    % 加密及构造ELEMENTs
    %%% ------------------------------ Refinement --------------------------
    Fracture{i}.EleFirst = countEle;
    for ii = 1 : nCrossPoint
        ni = Fracture{i}.eleNum(ii);
        if ni > 0.9
            theta = linspace(-pi,0,ni+1);
            % elem = (cos(theta)+1)*Fracture{i}.Li(ii)/2;
            elem = linspace(0,1,ni+1)*Fracture{i}.Li(ii);
            if ii < 1.1
                element = elem'*[Fracture{i}.cost,Fracture{i}.sint];
                element(:,1) = FracLoc(i,1) +element(:,1);
                element(:,2) = FracLoc(i,2) +element(:,2);
            else
                element = elem'*[Fracture{i}.cost,Fracture{i}.sint];
                element(:,1) = element(:,1) +Point(ii-1,1);
                element(:,2) = element(:,2) +Point(ii-1,2);
            end
            [len,~]= size(element);
            Fracture{i}.element(ii,1:len,1:2) = element;
            AllEle_global(countEle:countEle+len-2,1:4) = [element(1:len-1,:),element(2:len,:)];
            AllEle_global(countEle:countEle+len-2,5) = Fracture{i}.sint;
            AllEle_global(countEle:countEle+len-2,6) = Fracture{i}.cost;
            AllEle_global(countEle:countEle+len-2,7) = CalculateDis(element(2:len,:),element(1:len-1,:));
            AllEle_global(countEle:countEle+len-2,8) = (AllEle_global(countEle:countEle+len-2,1)+AllEle_global(countEle:countEle+len-2,3))/2;
            AllEle_global(countEle:countEle+len-2,9) = (AllEle_global(countEle:countEle+len-2,2)+AllEle_global(countEle:countEle+len-2,4))/2;
            countEle = countEle + len-1;
        end
    end
    if nCrossPoint > 0.1
        ii = ii+1;
        ni = Fracture{i}.eleNum(ii);
        theta = linspace(-pi,0,ni+1);
        % elem = (cos(theta)+1)*Fracture{i}.Li(ii)/2;
        elem = linspace(0,1,ni+1)*Fracture{i}.Li(ii);
        element = elem'*[Fracture{i}.cost,Fracture{i}.sint];
        element(:,1) = element(:,1) +Point(ii-1,1);
        element(:,2) = element(:,2) +Point(ii-1,2);
        [len,~]= size(element);
        AllEle_global(countEle:countEle+len-2,1:4) = [element(1:len-1,:),element(2:len,:)];
        AllEle_global(countEle:countEle+len-2,5) = Fracture{i}.sint;
        AllEle_global(countEle:countEle+len-2,6) = Fracture{i}.cost;
        AllEle_global(countEle:countEle+len-2,7) = CalculateDis(element(2:len,:),element(1:len-1,:));
        AllEle_global(countEle:countEle+len-2,8) = (AllEle_global(countEle:countEle+len-2,1)+AllEle_global(countEle:countEle+len-2,3))/2;
        AllEle_global(countEle:countEle+len-2,9) = (AllEle_global(countEle:countEle+len-2,2)+AllEle_global(countEle:countEle+len-2,4))/2;
        countEle = countEle + len-1;
        Fracture{i}.EleLast = countEle-1;
        Fracture{i}.element(ii,1:len,1:2) = element;
        
    else
        ni = Fracture{i}.eleNum;
        theta = linspace(-pi,0,ni+1);
        %refine
        %elem = (cos(theta)+1)/2*Fracture{i}.length;
        %not refine
        elem = linspace(0,1,ni+1)*Fracture{i}.length;
        
        element = elem'*[Fracture{i}.cost,Fracture{i}.sint];
        element(:,1) = FracLoc(i,1) +element(:,1);
        element(:,2) = FracLoc(i,2) +element(:,2);
        [len,~]= size(element);
        AllEle_global(countEle:countEle+len-2,1:4) = [element(1:len-1,:),element(2:len,:)];
        AllEle_global(countEle:countEle+len-2,5) = Fracture{i}.sint;
        AllEle_global(countEle:countEle+len-2,6) = Fracture{i}.cost;
        AllEle_global(countEle:countEle+len-2,7) = CalculateDis(element(2:len,:),element(1:len-1,:));
        AllEle_global(countEle:countEle+len-2,8) = (AllEle_global(countEle:countEle+len-2,1)+AllEle_global(countEle:countEle+len-2,3))/2;
        AllEle_global(countEle:countEle+len-2,9) = (AllEle_global(countEle:countEle+len-2,2)+AllEle_global(countEle:countEle+len-2,4))/2;
        countEle = countEle + len-1;
        Fracture{i}.element(1:len,1:2) = element;
        Fracture{i}.EleLast = countEle-1;
    end
    endn = countEle-1;
    %     CM_marker(endn) = 2;
    %     CM_marker(endn-1) = 1;
    %     CM_marker(endn-2) = 0;
    type = Fracture{i}.Type;
    AllEle_global(startn:endn,10) = type;
end
countEle = countEle-1;
if countEle~=nAllEle_global
    error('number of elements are not correct');
end
fprintf('The total number of elements are %d\n',nAllEle_global);
%Extra Connections
for i = 1 : nFracture
    nCrossPoint = Fracture{i}.nConnectF;
    for ii = 1 : nCrossPoint
        iindex = Fracture{i}.index(ii);
        NeighbFrac  = Fracture{i}.NeighbF(ii) ;
        nCrossPoint2 = Fracture{NeighbFrac}.nConnectF;
        if iindex == 38
            f = 1;
        end
        for jj = 1 :nCrossPoint2
            NeighbFrac2  = Fracture{NeighbFrac}.NeighbF(jj) ;
            if NeighbFrac2 == i
                jjindex = Fracture{NeighbFrac}.index(jj);
                % If there are no common points this connection is not
                % valid
                L = AllEle_global(iindex,:);
                dot = AllEle_global(jjindex,1:2);
                [sta1] = isDotIn(dot, L);
                dot = AllEle_global(jjindex,3:4);
                [sta2] = isDotIn(dot, L);
                if max(sta1,sta2) > 0.1
                    ConnList_global(iindex,2) = ConnList_global(iindex,2) + 1;
                    ConnList_global(iindex,ConnList_global(iindex,2)+2) = jjindex;
                end
                dot = AllEle_global(jjindex+1,1:2);
                [sta1] = isDotIn(dot, L);
                dot = AllEle_global(jjindex+1,3:4);
                [sta2] = isDotIn(dot, L);
                if max(sta1,sta2) > 0.1
                    ConnList_global(iindex,2) = ConnList_global(iindex,2) + 1;
                    ConnList_global(iindex,ConnList_global(iindex,2)+2) = jjindex+1;
                end
                
                L = AllEle_global(iindex+1,:);
                dot = AllEle_global(jjindex,1:2);
                [sta1] = isDotIn(dot, L);
                dot = AllEle_global(jjindex,3:4);
                [sta2] = isDotIn(dot, L);
                if max(sta1,sta2) > 0.1
                    ConnList_global(iindex+1,2) = ConnList_global(iindex+1,2) + 1;
                    ConnList_global(iindex+1,ConnList_global(iindex+1,2)+2) = jjindex;
                end
                dot = AllEle_global(jjindex+1,1:2);
                [sta1] = isDotIn(dot, L);
                dot = AllEle_global(jjindex+1,3:4);
                [sta2] = isDotIn(dot, L);
                if max(sta1,sta2) > 0.1
                    ConnList_global(iindex+1,2) = ConnList_global(iindex+1,2) + 1;
                    ConnList_global(iindex+1,ConnList_global(iindex+1,2)+2) = jjindex+1;
                end
            end
        end
    end
end
for iindex = 1 : nAllEle_global
    if AllEle_global(iindex,10) > 6 && ConnList_global(iindex,2) > 2
        AllEle_global(iindex,10)  = 1;
    end
    if AllEle_global(iindex,7) < 1e-3
        ConnList_global(iindex,2) = 0;
    end
    if ConnList_global(iindex,2) > 2
        cons = zeros(5,1)+1;
       
        nci = ConnList_global(iindex,2);
        nc0 = nci;
        for in = 1 : nci
            nb  = ConnList_global(iindex,2+in);
            if nb > -0.1
                point1 = AllEle_global(nb,1:2);
                isin = isDotIn(point1,AllEle_global(iindex,1:4));
                if isin < 0.1
                    point1 = AllEle_global(nb,3:4);
                    isin = isDotIn(point1,AllEle_global(iindex,1:4));
                end
                if isin < 0.1
                    cons(in) = 0;
                    %not connect
                    nc0 = nc0 - 1;
                    ConnList_global(iindex,2+in) = -1;
                end
                if AllEle_global(nb,7) < 1e-2
                    cons(in) = 0;
                    nc0 = nc0 - 1;
                    ConnList_global(iindex,2+in) = -1;
                end
            end
        end
        conv = zeros(5,1);
        num = 0;
        for inn = 1 : 5
            if cons(inn) > 0.1
                if ConnList_global(iindex,2+inn) > 0.1
                num = num + 1;
                conv(num) = ConnList_global(iindex,2+inn);
                end
            end
        end
        ConnList_global(iindex,3:7) = conv;
                
        ConnList_global(iindex,2) = nc0;
        %Reconstruct
        
        if nc0 > 3
            f = 1;
        end
    end
end
% %可视化剖分的网格
%  for i = 1 : nFracture
%      plot([FracLoc(i,1) FracLoc(i,3)],[FracLoc(i,2) FracLoc(i,4)],'c','Linewidth',6);
%      hold on;
%      nc = Fracture{i}.nConnectF;
%      if nc < 0.1
%          for ii = 1 : Fracture{i}.eleNum
%              plot(Fracture{i}.element(:,1),Fracture{i}.element(:,2),'.')
%              hold on;
%          end
%      else
%          for ii = 1 : nc+1
%              plot(Fracture{i}.element(ii,:,1),Fracture{i}.element(ii,:,2),'.')
%              hold on;
%          end
%      end
%  end
Fractures = Fracture;
nFractures = nFracture;
end
function fid=ReadComment(fid)
pos0 = ftell(fid);
temp = fscanf(fid,'%s',1);
while strcmp(temp,'--')==1
    temp = fgetl(fid);
    pos0 = ftell(fid);
    temp = fscanf(fid,'%s',1);
end
fseek(fid,pos0,'bof');
end