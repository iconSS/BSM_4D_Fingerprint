clear all;
clc;

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

Pairwise = {};
nPairWise = 0;

for s = 1:size(Infor_List,1)
    nPairWise = nPairWise+1;
    Pairwise{nPairWise,1} = Infor_List{s,1};
    Pairwise{nPairWise,2} = Infor_List{s,2};
    Pairwise{nPairWise,3} = Infor_List{s,3};
    
    Similarity_Total = 0;
    Compare_Total = 0;
    Compare_Total_Effective = 0;
    
    for i = 1:size(NameList_Train,1)
        for j = 1:size(NameList_Train,1)
            if i<j
                Name1 = NameList_Train{i,1};
                Name2 = NameList_Train{j,1};
                
                %%%%% same Name; Same group; same mode
                if strcmp(Name1(1:2),Infor_List{s,1}) && strcmp(Name2(1:2),Infor_List{s,1})
                    if strcmp(Name1(4:5),Infor_List{s,2}) && strcmp(Name2(4:5),Infor_List{s,2})
                        if strcmp([Name1(14),Name1(10)],Infor_List{s,3}) && strcmp([Name2(14),Name2(10)],Infor_List{s,3})
                            MatchedValue = Sim_CRB(i,j);
                            
                            Compare_Total = Compare_Total+1;
                            if MatchedValue>0.1
                                Similarity_Total = Similarity_Total+MatchedValue;
                                Compare_Total_Effective = Compare_Total_Effective+1;
                            end
                        end
                    end
                end
            end
        end
    end
    %                 Pairwise{nPairWise,3} = countMatch;
    Pairwise{nPairWise,4} = Similarity_Total;
    Pairwise{nPairWise,5} = Compare_Total_Effective;
    Pairwise{nPairWise,6} = Compare_Total;
    Pairwise{nPairWise,7} = (Pairwise{nPairWise,6}/(Pairwise{nPairWise,4}))^2/100; %%%% reversed average 1/Similarity
    Pairwise{nPairWise,8} = log(Pairwise{nPairWise,7}); %% MecDistance
end


