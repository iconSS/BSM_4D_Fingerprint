clear all;
clc

%%%%%%% calculate distance between sample Names
load D4_Sim_Train.mat

for Step1 = 1:1 %%%%%%% Summarize all the names groups modes of datasets
nSampleName = 0;
SampleName = {};
for i = 1:size(NameList_Train,1)
    if nSampleName == 0
        nSampleName = nSampleName+1;
        SampleName{nSampleName,1} = NameList_Train{i,1}(1:2);
    else
        vot = 0;
        for u = 1:nSampleName            
            if strcmp(NameList_Train{i,1}(1:2),SampleName{u,1}(1:2))
                vot = 1;
                break;
            end
        end
        if vot == 0
            nSampleName = nSampleName+1;
            SampleName{nSampleName,1} = NameList_Train{i,1}(1:2);
        end
    end
end

SampleName_Group = SampleName;
for s = 1:nSampleName
    nGroupName = 0;    
    for i = 1:size(NameList_Train,1)
        if strcmp(NameList_Train{i,1}(1:2),SampleName_Group{s,1}(1:2))
            if nGroupName == 0
                nGroupName = nGroupName+1;
                SampleName_Group{s,nGroupName+1} = NameList_Train{i,1}(4:5);
            else
                vot = 0;
                for u = 1:nGroupName
                    if strcmp(NameList_Train{i,1}(4:5),SampleName_Group{s,u+1}(1:2))
                        vot = 1;
                        break;
                    end
                end
                if vot == 0
                    nGroupName = nGroupName+1;
                    SampleName_Group{s,nGroupName+1} = NameList_Train{i,1}(4:5);
                end
            end
        end
    end
end

SampleName_Mode = SampleName;
for s = 1:nSampleName    
    nModeName = 0;    
    for i = 1:size(NameList_Train,1)
        if strcmp(NameList_Train{i,1}(1:2),SampleName_Mode{s,1}(1:2))
            if nModeName == 0
                nModeName = nModeName+1;
                SampleName_Mode{s,nModeName+1} = [NameList_Train{i,1}(14),NameList_Train{i,1}(10)];
            else
                vot = 0;
                for u = 1:nModeName
                    if strcmp(NameList_Train{i,1}(14),SampleName_Mode{s,u+1}(1)) && strcmp(NameList_Train{i,1}(10),SampleName_Mode{s,u+1}(2))
                        vot = 1;
                        break;
                    end
                end
                if vot == 0
                    nModeName = nModeName+1;
                    SampleName_Mode{s,nModeName+1} = [NameList_Train{i,1}(14),NameList_Train{i,1}(10)];
                end
            end
        end
    end
end

count_Infor = 0;

for s = 1:size(SampleName,1)
    for g = 2:size(SampleName_Group,2)
        if ~isempty(SampleName_Group{s,g})
            for m = 2:size(SampleName_Mode,2)
                if ~isempty(SampleName_Mode{s,m})
                    count_Infor = count_Infor+1;
                    Infor_List{count_Infor,1} = SampleName{s,1};
                    Infor_List{count_Infor,2} = SampleName_Group{s,g};
                    Infor_List{count_Infor,3} = SampleName_Mode{s,m};
                end
            end
        end
    end
end
end


for Step2 = 1:1 %%%%% calculate threshold value of similarity
Inter_Sim = zeros(count_Infor,7);
count_Infor = 0;
for s = 1:size(SampleName,1)
    t00 = clock;
    
    for g = 2:size(SampleName_Group,2)
        if ~isempty(SampleName_Group{s,g})
            for m = 2:size(SampleName_Mode,2)
                if ~isempty(SampleName_Mode{s,m})
                    count_Infor = count_Infor+1;
                    for i = 1:size(NameList_Train,1)
                        for j = 1:size(NameList_Train,1)
                            Name1 = NameList_Train{i,1};
                            Name2 = NameList_Train{j,1};
                            %%% Same Name; Same Group; Same Ion Mode && Same Chrom Mode
                            if strcmp(Name1(1:2),SampleName{s,1}) && strcmp(Name1(4:5),SampleName_Group{s,g}) ...
                                    && strcmp([Name1(14),Name1(10)],SampleName_Mode{s,m}) 
                                if strcmp(Name1(1:2),Name2(1:2)) && strcmp(Name1(4:5),Name2(4:5)) ...
                                        && strcmp(Name1(10:12),Name2(10:12)) && strcmp(Name1(14:15),Name2(14:15)) ...
                                        
                                    Inter_Sim(count_Infor,1) =  Inter_Sim(count_Infor,1)+Sim_CRB(i,j);
                                    Inter_Sim(count_Infor,2) = Inter_Sim(count_Infor,2)+1;
                                else
                                    Inter_Sim(count_Infor,3) =  Inter_Sim(count_Infor,3)+Sim_CRB(i,j);
                                    Inter_Sim(count_Infor,4) = Inter_Sim(count_Infor,4)+1;
                                end
                            end
                        end
                    end
                    
                end
            end
        end
    end
    
    t01=clock-t00;
    timetake=t01(1,6)+t01(1,5)*60+t01(1,4)*3600;
    disp(['Sample No.', num2str(s), '/','total No. ', num2str(size(SampleName,1)), ' takes ', num2str(timetake), ' seconds ']);
    
end

count_Infor = 0;
for s = 1:size(SampleName,1)
    for g = 2:size(SampleName_Group,2)
        if ~isempty(SampleName_Group{s,g})
            for m = 2:size(SampleName_Mode,2)
                if ~isempty(SampleName_Mode{s,m})
                    count_Infor = count_Infor+1;
                    
                    Inter_Sim(count_Infor,5) = Inter_Sim(count_Infor,1)/Inter_Sim(count_Infor,2);
                    Inter_Sim(count_Infor,6) = Inter_Sim(count_Infor,3)/Inter_Sim(count_Infor,4);
                    Inter_Sim(count_Infor,7) = Inter_Sim(count_Infor,5)/Inter_Sim(count_Infor,6);
                end
            end
        end
    end
end
end


Pairwise = {};
nPairWise = 0;
for Step3 = 1:1 %%%% comparison between Samples
    for s1 = 1:size(SampleName,1)
        for s2 = 1:size(SampleName,1)
            if s1<s2
                nPairWise = nPairWise+1;
                Pairwise{nPairWise,1} = SampleName{s1,1};
                Pairwise{nPairWise,2} = SampleName{s2,1};
                
                countMatch = 0; %%%% Number of matched fingerprints
                
                for i = 1:size(NameList_Train,1)
                    for j = 1:size(NameList_Train,1)
                        Name1 = NameList_Train{i,1};
                        Name2 = NameList_Train{j,1};
                        
                        %%%%% different Name; same mode
                        if strcmp(Name1(1:2),SampleName{s1,1}) && strcmp(Name2(1:2),SampleName{s2,1})
                            if strcmp(Name1(10),Name2(10)) && strcmp(Name1(14),Name2(14)) 
                                Threshold_Value = zeros(1,100);
                                nTV = 0;
                                for p = 1:size(Infor_List,1)
                                    if strcmp(Name1(1:2),Infor_List{p,1}) || strcmp(Name2(1:2),Infor_List{p,1})
                                        if strcmp([Name1(14),Name1(10)],Infor_List{p,3}) || strcmp([Name2(14),Name2(10)],Infor_List{p,3})
                                            nTV = nTV+1;
                                            Threshold_Value(nTV,1) = Inter_Sim(p,5);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end










Inter_Sim = zeros(4,4);
%%%% defined the different spectra higher than mean  
TP = zeros(4,1);
TN = zeros(4,1);
FP = zeros(4,1);
FN = zeros(4,1);
for i = 1:size(NameList_Train,1)    
    for j = 1:size(NameList_Train,1)
        Name1 = NameList_Train{i,1};
        Name2 = NameList_Train{j,1};
        if strcmp(Name1(1:2),Name2(1:2)) && strcmp(Name1(4:5),Name2(4:5)) && strcmp(Name1(10:12),Name2(10:12))
            Inter_Sim(1,1) =  Inter_Sim(1,1)+Sim_CRB(i,j);
            Inter_Sim(1,2) = Inter_Sim(1,2)+1;
        else
            Inter_Sim(1,3) =  Inter_Sim(1,3)+Sim_CRB(i,j);
            Inter_Sim(1,4) = Inter_Sim(1,4)+1;
        end
    end
end

for i = 1:1
    Inter_Sim(i,5) = Inter_Sim(i,1)/Inter_Sim(i,2);
    Inter_Sim(i,6) = Inter_Sim(i,3)/Inter_Sim(i,4);
    Inter_Sim(i,7) = Inter_Sim(i,5)/Inter_Sim(i,6);
end

for i = 1:size(NameList_Train,1)    
    for j = 1:size(NameList_Train,1)
        Name1 = NameList_Train{i,1};
        Name2 = NameList_Train{j,1};
        if strcmp(Name1(1:2),Name2(1:2)) && strcmp(Name1(4:5),Name2(4:5)) && strcmp(Name1(10:12),Name2(10:12))
            if Sim_CRB(i,j)>= Inter_Sim(1,6)
                TP(1,1) = TP(1,1)+1;
            else
                FN(1,1) = FN(1,1)+1;
            end
        end
    end
end

for i = 1:size(NameList_Train,1)
    for j = 1:size(NameList_Train,1)
        Name1 = NameList_Train{i,1};
        Name2 = NameList_Train{j,1};
        if strcmp(Name1(1:2),Name2(1:2)) && strcmp(Name1(4:5),Name2(4:5)) && strcmp(Name1(10:12),Name2(10:12))
        else            
            if Sim_CRB(i,j)>= Inter_Sim(1,5)
                FP(1,1) = FP(1,1)+1;
            else
                TN(1,1) = TN(1,1)+1;
            end
        end        
    end
end

for i = 1:1
    Inter_Sim(i,8) = TP(i,1)/(TP(i,1)+FP(i,1)); %TPR
    Inter_Sim(i,9) = TP(i,1)/(TP(i,1)+FN(i,1)); % Recall
    Inter_Sim(i,10) = FP(i,1)/(TP(i,1)+FN(i,1)); %FDR
    Inter_Sim(i,11) = (TP(i,1)+TN(i,1))/(TP(i,1)+TN(i,1)+FP(i,1)+FN(i,1)); % accuracy
end