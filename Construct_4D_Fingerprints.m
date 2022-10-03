clear all
clc

Analysis_Time_Range = 27;

load 'Untargeted_Ion_List.mat'
SampleNo= 1;

Mats=dir(fullfile('.','*.mat'));
[nOfMats,~]=size(Mats);

fileList = []; nfileList = 0;
for iMat = 1:nOfMats
   if length(Mats(iMat).name) == 18 && ...
           strcmp(Mats(iMat).name(4),'_') && ...
           strcmp(Mats(iMat).name(6),'_') && ...
           strcmp(Mats(iMat).name(10),'_')&& ...
           strcmp(Mats(iMat).name(13),'_')
       nfileList = nfileList+1;
       fileList(nfileList,1) = nfileList;
       fileList(nfileList,2) = iMat;
   end
end



t1=clock;
for iMat= fileList(SampleNo,2)
    if  length(Mats(iMat).name)>2
        for ry=1:100
            if (strcmp(Mats(iMat).name(ry:ry+3),'.mat'))
                ry_P=ry;
                break;
            end
        end

        MS=load(Mats(iMat).name);
        TIC = MS.TIC;
        Data = MS.SpotSummary;
        ST = MS.ScanTime;
        
        if (strcmp(Mats(iMat).name(ry_P-8:ry_P-6),'POS'))
            CRB_List = CRB_SWATH_Ion_List_RP;
            DIAL_List = DIAL_SWATH_Ion_List_RP;
            DDA_List = DDA_SWATH_Ion_List_RP;
        end
        if (strcmp(Mats(iMat).name(ry_P-8:ry_P-6),'NEG'))
            CRB_List = CRB_SWATH_Ion_List_RN;
            DIAL_List = DIAL_SWATH_Ion_List_RN;
            DDA_List = DDA_SWATH_Ion_List_RN;
        end
        
        CRB_List = CRB_List(CRB_List(:,2)>0 & CRB_List(:,2)<Analysis_Time_Range,:);
        
        DDA_List_Dereplicated = DDA_List;        
        DDA_List_Dereplicated(:,1) = round(DDA_List_Dereplicated(:,1)*10^2)/10^2;
        DDA_List_Dereplicated(:,2) = round(DDA_List_Dereplicated(:,2)*10^1)/10^1;
        DDA_List_Dereplicated = unique(DDA_List_Dereplicated, 'rows');

        for i = 1:size(DDA_List_Dereplicated,1)
            Temp = DDA_List(abs(DDA_List(:,1)-DDA_List_Dereplicated(i,1))<=0.01 & ...
                abs(DDA_List(:,2)-DDA_List_Dereplicated(i,2))<=0.1,:);
            Temp = sortrows(Temp,1);
            Median_mz = Temp(fix(size(Temp,1)/2)+1,1);
            Temp1 = Temp(Temp(:,1) == Median_mz,:);
            Median_tR = Temp1(fix(size(Temp1,1)/2)+1,2);
            DDA_List_Dereplicated(i,1) = Median_mz; 
            DDA_List_Dereplicated(i,2) = Median_tR; 
        end
        DDA_List_Dereplicated = unique(DDA_List_Dereplicated, 'rows');
        DDA_List_Dereplicated = DDA_List_Dereplicated(DDA_List_Dereplicated(:,1)>0,:);
        DDA_List = DDA_List_Dereplicated;
        
        
        t2=clock-t1;
        runtime=t2(1,6)+t2(1,5)*60+t2(1,4)*3600;
        if runtime<0
            runtime=runtime+3600*24;
        end
        disp(['Processing List ',' = ', num2str(runtime),'seconds']);
        
        t1 = clock;
        CRB_EIC_Extract = Extracting_EIC_CRB(CRB_List,Data,ST);
        t2=clock-t1;
        runtime=t2(1,6)+t2(1,5)*60+t2(1,4)*3600;
        if runtime<0
            runtime=runtime+3600*24;
        end
        disp([Mats(iMat).name,'CRB Processing ',' = ', num2str(runtime),'seconds']);
        
        t1 = clock;
        DIAL_EIC_Extract = Extracting_EIC(DIAL_List,Data,ST);
        t2=clock-t1;
        runtime=t2(1,6)+t2(1,5)*60+t2(1,4)*3600;
        if runtime<0
            runtime=runtime+3600*24;
        end
        disp([Mats(iMat).name,'DIAL Processing ',' = ', num2str(runtime),'seconds']);
        
        t1=clock;
        DDA_EIC_Extract = Extracting_EIC(DDA_List,Data,ST);
        t2=clock-t1;
        runtime=t2(1,6)+t2(1,5)*60+t2(1,4)*3600;
        if runtime<0
            runtime=runtime+3600*24;
        end
        disp([Mats(iMat).name,'DDA Processing ',' = ', num2str(runtime),'seconds']);
        
        [~, outputFile0] = fileparts(Mats(iMat).name(1:ry_P));
        outputFilename = strcat('Finger_',outputFile0,'.mat');        
        save(outputFilename,'CRB_EIC_Extract','DIAL_EIC_Extract','DDA_EIC_Extract','TIC','-v7.3');
        
        t2=clock-t1;
        runtime=t2(1,6)+t2(1,5)*60+t2(1,4)*3600;
        if runtime<0
            runtime=runtime+3600*24;
        end
        disp([Mats(iMat).name,' Processing ',' = ', num2str(runtime),'seconds']);
    end
end



function CRB_SWATH_Ion_List = CRB_Peak_Listing(datapoolfilepath,Dataprefix)
CRB_SWATH_Ion_List= [];
Mats=dir(fullfile(datapoolfilepath,'*.mat'));
[nOfMats,~]=size(Mats);
for iMat=1:nOfMats
    if  length(Mats(iMat).name)>size(Dataprefix,2) && (strcmp(Mats(iMat).name(1:size(Dataprefix,2)),Dataprefix))
        for ry=1:100
            if (strcmp(Mats(iMat).name(ry:ry+3),'.mat'))
                ry_P=ry;
                break;
            end
        end
        MS=load(Mats(iMat).name);
        IonList = MS.RB_Dereplicate_WithMSMS;
        IonList_Temp = IonList(:,2); %%% mz 
        IonList_Temp(:,2) = IonList(:,11); %%% mz 
        IonList_Temp = IonList_Temp(IonList_Temp(:,2)>0,:);
        CRB_SWATH_Ion_List = [CRB_SWATH_Ion_List;IonList_Temp];
    end
end
CRB_SWATH_Ion_List = sortrows(CRB_SWATH_Ion_List,1);
end

function DIAL_SWATH_Ion_List = DIAL_Peak_Listing(datapoolfilepath,Dataprefix)
DIAL_SWATH_Ion_List= [];
Mats=dir(fullfile(datapoolfilepath,'*.mat'));
[nOfMats,~]=size(Mats);
for iMat=1:nOfMats
    if  length(Mats(iMat).name)>size(Dataprefix,2) && (strcmp(Mats(iMat).name(1:size(Dataprefix,2)),Dataprefix))
        for ry=1:100
            if (strcmp(Mats(iMat).name(ry:ry+3),'.mat'))
                ry_P=ry;
                break;
            end
        end
        MS=load(Mats(iMat).name);
        IonList = MS.Exp_MSMS;
        IonList_Temp = zeros(size(IonList,1),2);
        for i = 1:size(IonList_Temp,1)
            IonList_Temp(i,1) = IonList{i,1};
            IonList_Temp(i,2) = IonList{i,2};
        end
        DIAL_SWATH_Ion_List = [DIAL_SWATH_Ion_List;IonList_Temp];
    end
end
DIAL_SWATH_Ion_List = sortrows(DIAL_SWATH_Ion_List,1);
end

function DDA_SWATH_Ion_List = DDA_Peak_Listing(datapoolfilepath,Dataprefix)
DDA_SWATH_Ion_List= [];
Mats=dir(fullfile(datapoolfilepath,'*.mat'));
[nOfMats,~]=size(Mats);
for iMat=1:nOfMats
    if  length(Mats(iMat).name)>size(Dataprefix,2) && (strcmp(Mats(iMat).name(1:size(Dataprefix,2)),Dataprefix))
        for ry=1:100
            if (strcmp(Mats(iMat).name(ry:ry+3),'.mat'))
                ry_P=ry;
                break;
            end
        end
        MS=load(Mats(iMat).name);
        IonList = MS.Q1_Infor;
        IonList = [IonList,[1:size(IonList,1)]'];
        Precursor = round(IonList(:,1)*10^2)/10^2;
        Precursor(:,2) = round(IonList(:,2)*10^1)/10^1;
        Precursor = unique(Precursor,'rows');
        IonList_S = zeros(size(Precursor,1),3);
        for i = 1:size(Precursor,1)
            Temp_List = IonList(abs(round(IonList(:,1)*10^2)/10^2)-Precursor(i,1) == 0 & ...
                abs(round(IonList(:,2)*10^1)/10^1)-Precursor(i,2) == 0,:);
            Index_TL = find(Temp_List(:,3) == max(Temp_List(:,3)));
            IonList_S(i,1:2) = Temp_List(Index_TL(1,1),1:2);
            IonList_S(i,3) = Temp_List(Index_TL(1,1),4);
        end

        MSMSOriginal = MS.DDA_MSMS_Summary;
        
        IonList_Temp = zeros(size(IonList_S,1),2);
        IonList = IonList_S;
        nIon_withMSMS = 0;
        for i = 1:size(IonList_Temp,1)
            MSMS = MSMSOriginal{IonList(i,3),1};
            MSMS = Dynamic_Background_Substraction(MSMS);
            if ~isempty(MSMS)
                nIon_withMSMS= nIon_withMSMS+1;
                IonList_Temp(nIon_withMSMS,1) = IonList(i,1);
                IonList_Temp(nIon_withMSMS,2) = IonList(i,2);
            end
        end
        DDA_SWATH_Ion_List = [DDA_SWATH_Ion_List;IonList_Temp];
    end
end
DDA_SWATH_Ion_List = sortrows(DDA_SWATH_Ion_List,1);
end

function EIC_Extract = Extracting_EIC_CRB(ion_infor,Data,ST)
for i = 1:size(ion_infor,1)
    tR = ion_infor(i,2);
    Index = find(abs(ST-tR)<2);
    if ~isempty(Index)
        EIC_Trace = zeros(size(Index,2),2);
        for j = 1:size(Index,2)
            if Index(1,j)<=size(Data,3)
                Temp = Data(:,:,Index(1,j));
                Temp = Temp(abs(Temp(:,1)-ion_infor(i,1))<0.015,:);
                EIC_Trace(j,1) = ST(1,Index(1,j));
                if ~isempty(Temp)
                    EIC_Trace(j,2) = sum(Temp(:,2));
                end
            end
        end
        EIC_Extract{i,1} = ion_infor(i,1); %%% mz
        EIC_Extract{i,2} = ion_infor(i,2); %%% tR
        EIC_Extract{i,3} = EIC_Trace;
        EIC_Extract{i,4} = ion_infor(i,3); %%% whether tR G5 calculated?
    end
end
end

function EIC_Extract = Extracting_EIC(ion_infor,Data,ST)
for i = 1:size(ion_infor,1)
    tR = ion_infor(i,2);
    Index = find(abs(ST-tR)<2);
    if ~isempty(Index)
        EIC_Trace = zeros(size(Index,2),2);
        for j = 1:size(Index,2)
            if Index(1,j)<=size(Data,3)
                Temp = Data(:,:,Index(1,j));
                Temp = Temp(abs(Temp(:,1)-ion_infor(i,1))<0.015,:);
                EIC_Trace(j,1) = ST(1,Index(1,j));
                if ~isempty(Temp)
                    EIC_Trace(j,2) = sum(Temp(:,2));
                end
            end
        end
        EIC_Extract{i,1} = ion_infor(i,1);
        EIC_Extract{i,2} = ion_infor(i,2);
        EIC_Extract{i,3} = EIC_Trace;
    end
end
end

function Temp = Dynamic_Background_Substraction(Temp)
MSMS_Merged = round(Temp(:,1)*10^2)/10^2;
MSMS_Merged = unique(MSMS_Merged,'rows');
Temp_S = zeros(size(MSMS_Merged,1),2);
for p = 1:size(MSMS_Merged,1)
    Temp0 = Temp(abs(Temp(:,1)-MSMS_Merged(p,1))<=0.01,:);
    Index0 = find(Temp0(:,2) == max(Temp0(:,2)));
    Temp_S(p,1) = Temp0(Index0(1,1),1);
    Temp_S(p,2) = round(Temp0(Index0(1,1),2)*10^2)/10^2;
end

Temp_S = unique(Temp_S,'rows');
Temp = Temp_S;
while size(Temp,1)>0
    IntList = unique(Temp(:,2),'rows');
    IntList = sortrows(IntList,-1);
    vot = 0;
    for p = 1:size(IntList,1)
        Temp1 = find(Temp(:,2) == IntList(p,1));
        if size(Temp1,1)>5
            Temp = Temp(Temp(:,2) > IntList(p,1),:);
            vot =1;
            break;
        end
    end
    if vot == 0
        BaseLineSpots_Index = find(Temp(:,2)==min(Temp(:,2)));
        count3 = 0; %%%%% 3 consecutive baseline value trigger a baseline elimination
        for g = 1:size(BaseLineSpots_Index,1)-1
            if BaseLineSpots_Index(g+1,1)-BaseLineSpots_Index(g,1) == 1
                count3 = count3+1;
            else
                count3 = 0;
            end
            if count3 == 2
                Temp = Temp(Temp(:,2)>min(Temp(:,2)),:);
                break;
            end
        end
        if size(BaseLineSpots_Index,1) == 1 || g == size(BaseLineSpots_Index,1)-1
            break; %%%% break while
        end
    end
end
end