function [Analytics_Layers,Run_Time,Close_Contact,Roadmap] = fastfindinfector(infectorID,Layers,Start_Time,End_Time,Durationtime,Dispersion)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% FASTFIND Quickly find suspicious Closer_Contact functions in batch      %
% *Version：2020-05-14 13:09 V6.0                                         %
% Users can collect log files based on the base station or AP information %
% in their area, and use the fastfindinfector analysis function to get    %
% close contact data. The basic data needs the comparison relationship    %
% between the people in the area and the location and description data    %
% of the base station.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% *Coding：Jason and James*
global T_Filetime;
global T_UserID;
global User_Num;
global BaseStation_AP_Num;
global BaseStation_AP_Info;
URL_GroundData = 'http://************************';
URL_KeptLog = 'http://*************************';
URL_Duration = 'http://*************************';
tic;
FCAT_Path = pwd;
f = waitbar(0,'Looking for files, please wait...');
pause(0.5);
if isfile(fullfile(FCAT_Path,'data','GroundData.mat')) == 1
    load(fullfile(FCAT_Path,'data','GroundData.mat'));
else
    Outfilename_Basic = websave(fullfile(FCAT_Path,'data','GroundData.mat'),URL_GroundData);
end
[Xrow,indexUser] = ismember(infectorID,T_UserID.UserID);
if sum(Xrow(:)) == 0
    Analytics_Layers = 0;
    Run_Time = toc;
    Close_Contact = table();
    Roadmap = table();
    close(f);
    return;
end
[~,Start_Locb] = ismember(Start_Time,string(T_Filetime.hour));
[~,End_Locb] = ismember(End_Time,string(T_Filetime.hour));
Inquiry_Table = table();
Inquiry_Table.KName = unique(eraseBetween(T_Filetime.KName(Start_Locb:End_Locb),24,26),'stable');
Inquiry_Table.SName = unique(eraseBetween(T_Filetime.SName(Start_Locb:End_Locb),21,23),'stable');
Inquiry_Table.Date = unique(eraseBetween(Inquiry_Table.KName,1,14),'stable');
Time_lookforfile = toc;
waitbar(.1,f,['Time to find files:',num2str(round(Time_lookforfile,2)),'s']);
pause(0.5);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Under the working path of the client, use "Query Log" to find whether the file data to be accessed%
% exists. If it cannot be found in the log, it will initiate a request to load the data to the
% server and load the data.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
tic;
for ii = 1:height(Inquiry_Table)
    if isfile(fullfile(FCAT_Path,'data',strcat(Inquiry_Table.KName{ii},'.mat'))) == 1
        load(fullfile(FCAT_Path,'data',Inquiry_Table.KName{ii}));
    else
        URL_KeptLog = strcat(URL_KeptLog,Inquiry_Table.Date{ii});
        Outfilename_TK = websave(fullfile(FCAT_Path,'data',strcat(Inquiry_Table.KName{ii},'.mat')),URL_KeptLog);
    end
    if isfile(fullfile(FCAT_Path,'data',strcat(Inquiry_Table.SName{ii},'.mat'))) == 1
        load(fullfile(FCAT_Path,'data',Inquiry_Table.SName{ii}));
    else
        URL_Duration = strcat(URL_Duration,Inquiry_Table.Date{ii});
        Outfilename_SD = websave(fullfile(FCAT_Path,'data',strcat(Inquiry_Table.KName{ii},'.mat')),URL_Duration);
    end
end
Time_download = toc;
waitbar(.4,f,['Time load files: ',num2str(round(Time_download,2)),'s'],' Start Analyzing: ^_^');
pause(0.5);
M = 1;
Possible(1,:) = indexUser(1,:);

if nargin < 6
    Dispersion = 0;
end
Time_Layers_Sum = 0;
Possible_ID_Dura = sparse(User_Num,AP_Num);
while 1
    % According to how many levels are input, determine the depth of analysis
    tic;
    if M > Layers
        break;
    end
    CurrentNum = length(indexUser(1,:));
        S_Analysis = sparse(User_Num,BaseStation_AP_Num);
        % Analyze users one by one and output their close contacts,the number of cycles ii means all
        % sparse matrix.
        for ii = Start_Locb:End_Locb
            eval(['S1 = S_Duration_',T_Filetime.hour{ii},';']);
            if nonzeros(S1(indexUser(1,kk),:)) ~= 0
                % Calculate coexistence duration, and re-assign all IDs other than the current query
                % ID to be less than or equal to query ID duration under different APs to ensure
                % coexistence-duration.
                if Dispersion ~= 0
                    k = find(S1);
                    S1(k) = S1(k).*(S1(k) >= Dispersion);
                end
                % For BaseStation or APs with indexUser corresponding to non-zero duration, keep their values,
                % and the durations of the remaining BaseStation or APs are set to 0;
                SA = S1(indexUser(1,kk),:);
                [~,existAP] = find(SA);
                SA1 = S1(:,existAP);
                S1 = sparse(User_Num,BaseStation_AP_Num);
                S1(:,existAP) = SA1;
                for jj = 1:length(nonzeros(SA))
                    [~,Locb] = find(SA);
                    SB = S1(:,Locb(1,jj));
                    SB(SB >= SA(:,Locb(1,jj))) = SA(:,Locb(1,jj));
                    S1(:,Locb(1,jj)) = SB;
                end
                S_Analysis = S_Analysis + S1;
            end
        end
        if (CurrentNum == 1 && isempty(nonzeros(S_Analysis)) == 1) || isempty(nonzeros(nonzeros(S_Analysis)>=Durationtime*60))
            Analytics_Layers = 0;
            Run_Time = toc;
            Close_Contact = table();
            Roadmap = table();
            close(f);
            return
        else
            Duration = round(sum(S_Analysis,2)/60);
            Possible_ID_Dura(kk,M) = Duration(indexUser(1,kk),1);
            
            if nonzeros(S_Analysis) ~= 0
                [UserRow,~] = find(S_Analysis >= Durationtime*60);
                UserRow = unique(UserRow,'stable');
                [~,Locb] = ismember(indexUser(1,kk),UserRow);
                UserRow(Locb,:) = [];
                if kk == 1
                    Possible(M+1,:) = 0;
                    Possible(M+1,1:length(UserRow)) = UserRow;
                else
                    UserRow_Num_C = length(nonzeros(Possible(M+1,:)));
                    UserRow_Num = length(UserRow)+UserRow_Num_C;
                    Possible(M+1,1:UserRow_Num) = [Possible(M+1,1:UserRow_Num_C),UserRow'];
                    Possible(M+1,1:length(unique(Possible(M+1,:)))) = unique(Possible(M+1,:),'stable');
                    Possible(M+1,length(unique(Possible(M+1,:)))+1:UserRow_Num) = 0;
                end
            end
    indexUser = nonzeros(Possible(M+1,:));
    P = Possible(1:M,:);
    [lia,~] = ismember(indexUser,P);
    indexUser = unique(indexUser(~lia)','stable');
    if isempty(indexUser) == 1
        break;
    end
    Possible(M+1,:) = [];
    Possible(M+1,1:length(indexUser)) = indexUser(1,:);
    Time_Layers = toc; % 输出分析时间
    Proportion = round(4/Layers*M+5);
    text1 = strcat('"No.',num2str(M),'Layer use time： ',num2str(round(Time_Layers,2)),'s"');
    eval(['waitbar(.',num2str(Proportion),',f,',text1,');']);
    pause(0.5);
    Time_Layers_Sum = Time_Layers_Sum+Time_Layers;
    M = M+1;
end
Possible(M,:) = [];
Possible = Possible';
Analytics_Layers = M-1;
Close_Contact = table();
for ii = 1:M-1
    Layers_Num = ii*ones(length(nonzeros(Possible(:,ii))),1);
    Possible_Table_Temp = table();
    Possible_Table_Temp.UserIDX= nonzeros(Possible(:,ii));
    Possible_Table_Temp.Layer = Layers_Num;
    Close_Contact = [Close_Contact;Possible_Table_Temp];
end
Close_Contact.UserID= T_UserID.UserID(Close_Contact.UserIDX);
Close_Contact.Duration = nonzeros(Possible_ID_Dura);
Close_Contact.Unit1 = T_UserID.Unit1(Close_Contact.UserIDX);
Close_Contact.Identity = T_UserID.Identity(Close_Contact.UserIDX);
Close_Contact.Unit2 = T_UserID.Unit2(Close_Contact.UserIDX);

Time_Report = toc;
waitbar(1,f,['Time of Create Report： ',num2str(round(Time_Report,2)),'s']);
pause(0.5);
close(f);
Run_Time = num2str(round(Time_lookforfile+Time_download+Time_Layers_Sum+Time_Report,2));
end