%% Requirements:
    % 2 files: 
        % 1) sleep edf file, 
        % 2) mat file of Kubios output of ECG R peaks 

%% NOTE: 
    % update in each place with "update here" (use the search function to
    % find all of these places)
    % if there are fewer than 32 channels, you'll need to update edf2mat_samefss  
    
%% Parameters

clear all; close all;

segmin=3; % at least 3 minutes for a stable stage
hrbwin=20; % window size around HR burst
avbin=5;  % bin size for averaging EEG powers
m=1.25; % threshold for HR burst detection mean(RR)-m*std(RR)
oc=0; %1 / 0 outlier control on eeg powers
swa=[0.1 4]; alpha=[8 13]; theta=[4 7]; sigma=[12 16]; % EEG bands
lf=[0.04 0.15]; hf=[0.15 0.4]; % LF, HF bands
swa1=[0.5 1];
swa2=[1 4];
RR_upperlim = 2;
RR_lowerlim = 0.55; 

% update here: with paths of the three file types
edfFolder = '';
hrvFolder = '';

% update here: with path of the desired output directory of the files
%outputFolder = '';

%% Begin Analysis Code

% grabbing all edf files (extracting values)
edfFileList = dir(fullfile(edfFolder, '*.edf'));  
hrvFileList = dir(fullfile(hrvFolder, '*.mat'));  

for edfFile = 1:numel(edfFileList)
    subVisit = edfFileList(edfFile).name; %grab beginning of file
    edfPathName =  [edfFileList(edfFile).folder filesep];
    edfFileName =  [edfFileList(edfFile).name];
    edfFullPath = [edfFileList(edfFile).folder filesep edfFileList(edfFile).name];
    % find corresponding kubios (hrv) 
    fileBegin = subVisit(1:7); % update here: update based on file naming convention (you'll want to ...
    % have the output of this line be the name to search in other folders
    % for the other corresponding files). 
    hrvIndex = find(contains({hrvFileList.name}, fileBegin)); % file index of HRV that we want in the list
    hrvFile = hrvFileList(hrvIndex);

    if isempty(hrvFile) % check HRV files are empty
        disp(['Skipping empty file: ', fileBegin]);
        continue; % skip to next file
    else
        load(strcat(hrvFile.folder, filesep, hrvFile.name)); % load HRV file
        disp(['RRs loaded for: ', fileBegin])
        [X,Channels,fs]=edf2mat_samefss(edfFullPath); % loading the EDF file
        disp(['EDFs loaded for: ', fileBegin])

        % EEG channel selection
        toChoose = {'F3' 'F4' 'C3' 'C4' 'O1' 'O2' 'Cz'}; % update here: if you want other channels
        cn1=ismember(Channels,toChoose); 
        Selection=find(cn1==1);
        XEEG=X(Selection,:);
        Xothers=X;Xothers(Selection,:)=[]; clear X;

        % filtering 
        disp('filtering...');
        [ba2,aa2]=butter(4,swa(1)/(fs/2),'high');
        [ba1,aa1]=butter(4,swa(2)/(fs/2),'low');
        X_delta=filtfilt(ba2,aa2,XEEG'); X_delta=filtfilt(ba1,aa1,X_delta); 
        X_deltaPWR=abs(hilbert(X_delta)); X_delta=X_delta'; X_deltaPWR=X_deltaPWR';
        
        [ba2,aa2]=butter(4,swa1(1)/(fs/2),'high');
        [ba1,aa1]=butter(4,swa1(2)/(fs/2),'low');
        X_delta1=filtfilt(ba2,aa2,XEEG'); X_delta1=filtfilt(ba1,aa1,X_delta1); 
        X_delta1PWR=abs(hilbert(X_delta1)); X_delta1=X_delta1'; X_delta1PWR=X_delta1PWR';
        
        [ba2,aa2]=butter(4,swa2(1)/(fs/2),'high');
        [ba1,aa1]=butter(4,swa2(2)/(fs/2),'low');
        X_delta2=filtfilt(ba2,aa2,XEEG'); X_delta2=filtfilt(ba1,aa1,X_delta2); 
        X_delta2PWR=abs(hilbert(X_delta2)); X_delta2=X_delta2'; X_delta2PWR=X_delta2PWR';
        
        [ba4,aa4]=butter(4,alpha(1)/(fs/2),'high');
        [ba3,aa3]=butter(4,alpha(2)/(fs/2),'low'); 
        X_alpha=filtfilt(ba4,aa4,XEEG'); X_alpha=filtfilt(ba3,aa3,X_alpha); 
        X_alphaPWR=abs(hilbert(X_alpha)); X_alpha=X_alpha'; X_alphaPWR=X_alphaPWR';
        
        [ba6,aa6]=butter(4,theta(1)/(fs/2),'high');
        [ba5,aa5]=butter(4,theta(2)/(fs/2),'low');
        X_theta=filtfilt(ba6,aa6,XEEG'); X_theta=filtfilt(ba5,aa5,X_theta);
        X_thetaPWR=abs(hilbert(X_theta)); X_theta=X_theta'; X_thetaPWR=X_thetaPWR';
        
        [ba8,aa8]=butter(4,sigma(1)/(fs/2),'high');
        [ba7,aa7]=butter(4,sigma(2)/(fs/2),'low'); 
        X_sigma=filtfilt(ba8,aa8,XEEG'); X_sigma=filtfilt(ba7,aa7,X_sigma); 
        X_sigmaPWR=abs(hilbert(X_sigma)); X_sigma=X_sigma'; X_sigmaPWR=X_sigmaPWR';
        %
        % eval(['x' Channels{1, 6}  '= X(6,:);']);
        
        if ~exist('RES')
            RES=Res;
        end
        clear Res
        fgr=fieldnames(RES.CNT.rate);
        fs2=eval(['RES.CNT.rate' '.' fgr{1,1}]);
        RR_tot_ind=round((RES.HRV.Data.T_RR-RES.CNT.Offset)*fs2);
        RR_tot=(RR_tot_ind(2:end)-RR_tot_ind(1:end-1))./fs2;
        RR_tot_ind(1)=[];
        RR_tot_time=RR_tot_ind./fs2;
        RRts=spline(RR_tot_ind./fs2,RR_tot,1/fs2:1/fs2:RR_tot_ind(end)/fs2);
        % RR filter
        [blf2,alf2]=butter(4,lf(2)/(fs/2),'low');
        [blf1,alf1]=butter(4,lf(1)/(fs/2),'high');
        RRlf=filtfilt(blf2,alf2,RRts);  RRlf=filtfilt(blf1,alf1,RRlf);
        
        [bhf2,ahf2]=butter(4,hf(1)/(fs/2),'high');
        [bhf1,ahf1]=butter(4,hf(2)/(fs/2),'low');
        RRhf=filtfilt(bhf2,ahf2,RRts);  RRhf=filtfilt(bhf1,ahf1,RRhf);
        hfpw=abs(hilbert(RRhf));
        %
        if fs2~=fs
            disp(['Warning! is sampling rate ' fs ' or ' fs2 '?']);
        end
        
        % make fake marker file (this creates a sleep scoring file with all
        % 7s, so that the rest will run properly). 
        num_epochs = floor(length(XEEG)/fs/30); % to get # of mrkr
        mrkr = 7 * ones(num_epochs, 1);
        mrkr(mrkr==-1)=7;
        t=find((mrkr(2:end)-mrkr(1:end-1))~=0);
        smp=[0;t*30*fs;size(mrkr,1)*30*fs]';
        bnd=zeros(length(smp)-1,2);
        for i=1:length(smp)-1
            bnd(i,:)=[smp(i)+1 smp(i+1)]; % beginning and end of each bout
        end
        Stage=mrkr([1;t+1]); % bout sleep stage
        % assignin('base','bnd',bnd);
        
        duration=(bnd(:,2)-bnd(:,1)+1)/(fs*60); 
        % bouts longer than segmin at each stage

        in7=find(Stage==7 & duration>segmin & bnd(:,2)<length(RRts));  % fake score number (can ignore)
        % assignin('base','RRts',RRts);
        
        disp('bouts extracted');
        % HRB analysis
        disp('HR burst analysis...')
        [sbj_hrb_ind7,sbj_hrb7,sbj_hrbHFPW7]=myhrbwindowedd(in7,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim); sbj_hrb_HFpwr_bins7=myhrbbinnedd_EEG(sbj_hrbHFPW7,avbin,fs,oc);
        sbj_hrb_RR_bins7=myhrbbinnedd_EEG(sbj_hrb7,avbin,fs,oc);
        sbj_hrb_EEG7=myhrbwindowedd_EEG(XEEG,hrbwin,fs,oc,sbj_hrb_ind7); 
        sbj_hrb_deltaPWR7=myhrbwindowedd_EEG(X_deltaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_deltapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_deltaPWR7,avbin,fs,oc);
        sbj_hrb_alphaPWR7=myhrbwindowedd_EEG(X_alphaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_alphapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_alphaPWR7,avbin,fs,oc);
        sbj_hrb_thetaPWR7=myhrbwindowedd_EEG(X_thetaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_thetapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_thetaPWR7,avbin,fs,oc);
        sbj_hrb_sigmaPWR7=myhrbwindowedd_EEG(X_sigmaPWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_sigmapwr_bins7=myhrbbinnedd_EEG(sbj_hrb_sigmaPWR7,avbin,fs,oc);
        sbj_hrb_delta1PWR7=myhrbwindowedd_EEG(X_delta1PWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_delta1pwr_bins7=myhrbbinnedd_EEG(sbj_hrb_delta1PWR7,avbin,fs,oc);
        sbj_hrb_delta2PWR7=myhrbwindowedd_EEG(X_delta2PWR,hrbwin,fs,oc,sbj_hrb_ind7); sbj_hrb_delta2pwr_bins7=myhrbbinnedd_EEG(sbj_hrb_delta2PWR7,avbin,fs,oc);
        sbj_hrb_ind7_quartile = zeros(size(sbj_hrb_ind7));

        
        %%%%
        disp('calculating bin data...');
        eval('ace.Channels= Channels(1, Selection);');
        for i=1:length(Selection)
            eval(['ace.binDelta_' Channels{1, Selection(i)}  '= sbj_hrb_deltapwr_bins7{i,1};']);
            eval(['ace.binSlowDelta_' Channels{1, Selection(i)}  '= sbj_hrb_delta1pwr_bins7{i,1};']);
            eval(['ace.binFastDelta_' Channels{1, Selection(i)}  '= sbj_hrb_delta2pwr_bins7{i,1};']);
            eval(['ace.binAlpha_' Channels{1, Selection(i)}  '= sbj_hrb_alphapwr_bins7{i,1};']);
            eval(['ace.binSigma_' Channels{1, Selection(i)}  '= sbj_hrb_sigmapwr_bins7{i,1};']);
            eval(['ace.binTheta_' Channels{1, Selection(i)}  '= sbj_hrb_thetapwr_bins7{i,1};']);
        end
        eval(['ace.binHF_nostage'  '= sbj_hrb_HFpwr_bins7;']);
        eval(['ace.binRR_nostage'  '= sbj_hrb_RR_bins7;']);
        %%%%sleep
        disp('averaging...');%(find(sbj_hrb_ind2_quartile ==1),:)
        for i=1:length(Selection)
            eval(['ace.avDelta_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_deltapwr_bins7{i,1}, 1:size(sbj_hrb_deltapwr_bins7{i,1},1));']);
            eval(['ace.avSlowDelta_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta1pwr_bins7{i,1}, 1:size(sbj_hrb_delta1pwr_bins7{i,1},1));']);
            eval(['ace.avFastDelta_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_delta2pwr_bins7{i,1}, 1:size(sbj_hrb_delta2pwr_bins7{i,1},1));']);
            eval(['ace.avAlpha_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_alphapwr_bins7{i,1}, 1:size(sbj_hrb_alphapwr_bins7{i,1},1));']);
            eval(['ace.avSigma_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_sigmapwr_bins7{i,1}, 1:size(sbj_hrb_sigmapwr_bins7{i,1},1));']);
            eval(['ace.avTheta_' Channels{1, Selection(i)}  '= mean_rj_olyr(sbj_hrb_thetapwr_bins7{i,1}, 1:size(sbj_hrb_thetapwr_bins7{i,1},1));']);
            
        end

        eval(['ace.avHF_nostage'  '= mean_rj_olyr(sbj_hrb_HFpwr_bins7, 1:size(sbj_hrb_HFpwr_bins7,1));']);      
        eval(['ace.RRtimeseries'  '= RRts;']);
        eval(['ace.HRB_idx'  '= sbj_hrb_ind7;']);
        eval(['ace.HRB_quartile'  '= sbj_hrb_ind7_quartile;']);
        eval(['ace.HRB_density'  '= length(sbj_hrb_ind7)/sum(duration(in7));']);
        eval(['ace.HRB_nostage'  '= sbj_hrb7;']);
        eval(['ace.HRB_EEG_allCh'   '= sbj_hrb_EEG7;']);
        eval(['ace.HRB_DeltaHilbAmp_allCh'   '= sbj_hrb_deltaPWR7;']);     
        eval(['ace.HRB_SlowDeltaHilbAmp_allCh'   '= sbj_hrb_delta1PWR7;']);
        eval(['ace.HRB_FastDeltaHilbAmp_allCh'   '= sbj_hrb_delta2PWR7;']);
        eval(['ace.HRB_SigmaHilbAmp_allCh'   '= sbj_hrb_sigmaPWR7;']);
        eval(['ace.HRB_ThetaHilbAmp_allCh'   '= sbj_hrb_thetaPWR7;']);
        eval(['ace.HRB_AlphaHilbAmp_allCh'   '= sbj_hrb_alphaPWR7;']);
        
        disp('non-Ace analysis...')
        inc7=myblactivityy(in7,bnd,sbj_hrb_ind7,hrbwin,fs);
        
        for i=1:length(Selection)
            eval(['ace.blDelta_' Channels{1, Selection(i)}  '= mean(X_deltaPWR(i,inc7));']);
            eval(['ace.blSlowDelta_' Channels{1, Selection(i)}  '= mean(X_delta1PWR(i,inc7));']);
            eval(['ace.blFastDelta_' Channels{1, Selection(i)}  '= mean(X_delta2PWR(i,inc7));']);
            eval(['ace.blAlpha_' Channels{1, Selection(i)}  '= mean(X_alphaPWR(i,inc7));']);
            eval(['ace.blSigma_' Channels{1, Selection(i)}  '= mean(X_sigmaPWR(i,inc7));']);
            eval(['ace.blTheta_' Channels{1, Selection(i)}  '= mean(X_thetaPWR(i,inc7));']);           
        end

            eval(['ace.blHF_' Channels{1, Selection(i)}  '= mean(hfpw(inc7));']);
            eval(['ace.blRR_' Channels{1, Selection(i)}  '= mean(RRts(inc7));']);
        
            ace.fs=fs;
            ace.filename=edfFileName;
            % update here: with desired output filename
            AceFile = [outputFolder 'ace_' edfFileName(1:end-4)];
            save(AceFile ,'ace', '-v7.3')
    end    
end

%% Functions

% function to read in the EDF file format (update here: if there are fewer
% than 32 channels)
function [X,Channels,fs]=edf2mat_samefss(fileloc)
    % [FileName,PathName] = uigetfile('/bazhlab/naji/home/EDFs_ACH_500Hz/*.edf','Select the edf data file');
    % load([PathName FileName]);
    fid=fopen(fileloc);% in format of [PathName FileName]
    a=fread(fid,236,'*char');
    ndr=fread(fid,8,'*char');
    ndr=str2double(ndr'); %number of data records in sec
    a=fread(fid,8,'*char');
    drdur=str2double(a'); %duration of each data record in sec
    ns=fread(fid,4,'*char'); ns=ns'; ns=str2double(ns);% number of signal channels
    Channels=cell(1,ns);
    for i=1:ns
        C=fread(fid,16,'*char');C=C';
        Channels{i}=C(find(isspace(C)==0));
    end

    Channels = {Channels{1:32}}; % update here: Remove the last channels are not being saved ('Abdo'	'Light'	'PositionSen'	'SpO2'	'Pulse'	'OxStatus'	'Pleth'	'SCL'	'DerivedHR'	'EDFAnnotations')
    
    fread(fid,ns*80,'*char'); % channel transducer type can be extracted
    fread(fid,ns*8,'*char'); %channel physical dimension can be extracted
    phmn=zeros(1,ns); phmx=phmn;dmn=phmn;dmx=dmn;
    for i=1:ns
        pm=fread(fid,8,'*char');pm=pm';
        phmn(i)=str2double(pm);
    end                         %phys min
    for i=1:ns
        pm=fread(fid,8,'*char'); pm=pm';
        phmx(i)=str2double(pm);%phys max
    end
    for i=1:ns
        dm=fread(fid,8,'*char');dm=dm'; 
        dmn(i)=str2double(dm);
    end                         %dig min
    for i=1:ns
        dx=fread(fid,8,'*char'); dx=dx';
        dmx(i)=str2double(dx);
    end                         %dig max
    scalefac=(phmx-phmn)./(dmx-dmn);
    dc=phmx-scalefac.*dmx;
    
    fread(fid,ns*80,'*char'); % prefilters
    nr=zeros(1,ns);
    for i=1:ns
        nrc=fread(fid,8,'*char'); nrc=nrc'; 
        nr(i)=str2double(nrc); %number of samples in each data record
    end
    % Exclude the last channel from further processing
    nr = nr(1:32); %  update here: Remove the last channels are not being saved ('Abdo'	'Light'	'PositionSen'	'SpO2'	'Pulse'	'OxStatus'	'Pleth'	'SCL'	'DerivedHR'	'EDFAnnotations')
    ns = 32; % update here: with last channel index that you would like to save.
    if all(nr == nr(1))
        fs = nr(1); % Sampling frequency
    else
        sprintf('Data cant be stored in a single matrix')
    end
    fread(fid,ns*32,'*char');
    ch_fs=nr/drdur;
    if mean(nr)==nr(1) && mean(ch_fs)==ch_fs(1)
    X=zeros(ns,nr(1)*ndr);
    % for i=1:ns
    %     X{i,1}=zeros(1,nr(i)*ndr);
    % end
    fs=ch_fs(1);
    end
    spins={'\\','|','/','-'};
        reverseStr = 'Reading EDF file ... ';
    for i=1:ndr
        for j=1:ns
            s=fread(fid,nr(j),'int16').*scalefac(j)+dc(j);s=s';
            X(j,(i-1)*nr(j)+1:i*nr(j))=s;
        end
        
    %     si=mod(i,4); 
    %     if si==0
    %         si=4;
    %     end
    %     msg = (spins{si});
    %     fprintf([reverseStr, msg]);
    %     reverseStr = repmat(sprintf('\b'), 1, 1);
        
    end
    fprintf('\n');
    fclose(fid);
end

function [sbj_hrb_ind,sbj_hrb,sbj_hrbHFPW]=myhrbwindowedd(in,bnd,fs,hrbwin,RRts,hfpw,m,RR_lowerlim, RR_upperlim)
sbj_hrb_ind=[];
    sbj_hrb=[];
    sbj_hrbHFPW=[];
if ~isempty(in)
    
    for i=1:length(in)
        pmin=myHRB_finder(RRts(bnd(in(i),1):bnd(in(i),2)),fs,m, RR_lowerlim, RR_upperlim);
        for j=1:length(pmin)
            sbj_hrb_ind=[sbj_hrb_ind;pmin(j)+bnd(in(i),1)-1];
        end
    end
    if ~isempty(sbj_hrb_ind)
        rj=find(sbj_hrb_ind/fs<(hrbwin/2) | (length(RRts)-sbj_hrb_ind)/fs<(hrbwin/2));
        sbj_hrb_ind(rj)=[];
    end
    for i=1:length(sbj_hrb_ind)
        sbj_hrb=[sbj_hrb;RRts(sbj_hrb_ind(i)-floor(hrbwin/2*fs)+1:sbj_hrb_ind(i)+floor(hrbwin/2*fs))];
        sbj_hrbHFPW=[sbj_hrbHFPW;hfpw(sbj_hrb_ind(i)-floor(hrbwin/2*fs)+1:sbj_hrb_ind(i)+floor(hrbwin/2*fs))];
    end
end
end

function inc=myblactivityy(in,bnd,hrbind,hrbwin,fs)
inc=[];
    for i=1:length(in)
        inc=[inc bnd(in(i),1):bnd(in(i),2)];
    end
    rj=[];
    for i=1:length(hrbind)
        rj=[rj hrbind(i)-floor(hrbwin/2*fs)+1:hrbind(i)+floor(hrbwin/2*fs)];
    end
    inc(find(ismember(inc,rj)==1))=[];
end

function hrb_EEG=myhrbwindowedd_EEG(X,hrbwin,fs,oc,sbj_hrb_ind)
hrb_EEG=cell(size(X,1),1);
for j=1:size(X,1)
    hrb_EEG{j,1}=zeros(length(sbj_hrb_ind),hrbwin*fs);
    
    for i=1:length(sbj_hrb_ind)
        hrb_EEG{j,1}(i,:)=X(j,sbj_hrb_ind(i)-floor(hrbwin/2*fs)+1:sbj_hrb_ind(i)+floor(hrbwin/2*fs));
    end
    if oc==1
    mm=max(hrb_EEG{j,1}');
    rj=find(mm>(mean(mm)+3*std(mm)));
    hrb_EEG{j,1}(rj,:)=[];
    end
end
    
end

function hrb_EEG_bins=myhrbbinnedd_EEG(x,avbin,fs,oc)
if iscell(x)
hrb_EEG_bins=cell(size(x,1),1);
for j=1:size(x,1)
    hrb_EEG_bins{j,1}=zeros(size(x{j,1},1),size(x{j,1},2)/(fs*avbin));
    
    for i=1:size(x{j,1},1)
        for ii=1:size(x{j,1},2)/(fs*avbin)
        hrb_EEG_bins{j,1}(i,ii)=mean(x{j,1}(i,(ii-1)*avbin*fs+1:ii*avbin*fs));
        end
    end
    if oc==1
        mm=max(hrb_EEG_bins{j,1}');
        rj=find(mm>(mean(mm)+3*std(mm)));
        hrb_EEG_bins{j,1}(rj,:)=[];
    end
end
else
    hrb_EEG_bins=[];
    for i=1:size(x,1)
        for ii=1:size(x,2)/(fs*avbin)
            hrb_EEG_bins(i,ii)=mean(x(i,(ii-1)*avbin*fs+1:ii*avbin*fs));
        end
    end
    if oc==1
        mm=max(hrb_EEG_bins');
        rj=find(mm>(mean(mm)+3*std(mm)));
        hrb_EEG_bins(rj,:)=[];
    end
end
 
end

function mn = mean_rj_olyr(x,idx)
if length(idx) == 0
    mn = NaN;
else
    if size(x(idx,:),1) == 0
        mn = NaN;

    elseif size(x(idx,:),1) == 1
        mn = x(idx,:);
    else
        x = x(idx,:);
        mm=max(x');
        rjc=find(mm>(mean(mm)+3*std(mm)));
        x(rjc,:)=[];
        mn = mean(x, 1);
    end
end
end

function pmin=myHRB_finder(samples,fs,m, RR_lowerlim, RR_upperlim)
% m=1.25;
allmin=find((samples(1:end-2)-samples(2:end-1))>=0 & (samples(2:end-1)-samples(3:end))<=0);
th=mean(samples(allmin))-m*std(samples(allmin)); %threshold
pmin=find((samples(1:end-2)-samples(2:end-1))>=0 & (samples(2:end-1)-samples(3:end))<=0 & samples(1:end-2)<=th);
jj=1; pmnsp=[]; rjcts=[]; cmpmn=[];

            for ii=1:length(pmin)-1
                if (pmin(ii+1)-pmin(ii))>10*fs
                    jj=jj+1;
                    pmnsp=[];
                end
                if (pmin(ii+1)-pmin(ii))<10*fs
                    pmnsp=[pmnsp pmin(ii) pmin(ii+1)];
                    cmpmn{jj,1}=pmnsp; 
                end
            end
            [rcmp,~]=size(cmpmn);
            for ii=1:rcmp
                if ~isempty(cmpmn{ii,1})
                    tempv=cmpmn{ii,1};
                    [~,kkn]=min(samples(tempv));
                    tempv(kkn)=[];
                    rjcts=[rjcts tempv];
                end
            end
            rjsmp=[];
            for ii=1:length(rjcts)
                rjsmp=[rjsmp find(pmin== rjcts(ii))];
            end
            pmin(rjsmp)=[];
%             falseR=find(samples(pmin)<0.55); pmin(falseR)=[];
            
            rj = [];
            for ii=1:length(pmin)
                win = samples(max(1,pmin - 10*fs): min(pmin+10*fs , length(samples)));
                if min(win) < RR_lowerlim || max(win)>RR_upperlim
                    rj=[rj;ii];
                end
            end
            pmin(rj)=[];
end
