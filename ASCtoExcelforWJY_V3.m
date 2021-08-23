clear
ASCfile=dir(fullfile('*.ASC'));        %读入目录下文件信息存储为结构体形式
ASCstr=struct2cell(ASCfile);  %将格式转为cell形式
ASCname=ASCstr(1,:);        %取出其中文件名单
[~,nASCfile]=size(ASCname);   %计算文件个数
OutputName=strcat('auto',ASCname,'.xlsx');

for iASCfile=1:1:nASCfile
    if strfind(ASCname{iASCfile},'.ASC')    %如果是xls文件格式 注意括号要使用cell的括号
        [ASC{iASCfile}]=importdata(ASCname{iASCfile}); 
    end
    ASCcell=[];
    ASCnumtocell=[];
    ASCstrtocell=[];
    ASCraw=[];
    ASCstr2double=[];
    ASCcell2mat=[];
    ASCcell=struct2cell(ASC{1,iASCfile});
    ASCnumtocell=num2cell(ASCcell{1,1});
    ASCstrtocell=ASCcell{2,1};
    [numrow,numcol]=size(ASCnumtocell);
    [strrow,strcol]=size(ASCstrtocell);
    rawrow=min(numrow,strrow);
    ASCraw=[ASCstrtocell(1:rawrow,:),ASCnumtocell(1:rawrow,:)]; 
    NaNx=cellfun(@isnan, ASCraw ,'UniformOutput',false);
    [NaNrow,~]=find(cellfun(@(x) isequal(1,x),NaNx(:,18)));
    ASCraw(NaNrow,:)=[];
%     test=ASCraw(:,1:2);
%     test1=str2double(test);
%     test2=str2double(ASCraw);
ASCstr2double=str2double(ASCraw(:,1:strcol));
ASCcell2mat=cell2mat(ASCraw(:,strcol+1:18));
ASCdelNaN{1,iASCfile}=[ASCstr2double,ASCcell2mat];
end

[~,nASCdelNaN]=size(ASCdelNaN);
for ifile=1:1:nASCdelNaN
RawTable=ASCdelNaN{1,ifile};
group=RawTable(:,9);
ind=find(arrayfun(@(x) isequal(0,x),group));
%delete RiseTime>3
%if ind(end)~=length(group)
for iind=1:1:length(ind)
RawTable(ind(iind)-(iind-1),:)=[];
end
%else
% for iind=1:1:length(ind)-1
% RawTable(ind(iind),:)=[]; 
% end

%Time modification
%RawTable(:,2)=1000*RawTable(:,2);
RawTable(:,2)=RawTable(:,2);
%Timesort
Table{1,ifile}=sortrows(RawTable,2);

%correlation analysis
RiseTime=Table{1,ifile}(:,4);
Amplitude=Table{1,ifile}(:,3);
K{1,ifile}=polyfit(RiseTime,Amplitude,1);
corfficientcorrelation{1,ifile}=K{1,ifile}(1,1);
%Frequency
[mTable,nTable]=size(Table{1,ifile});
T=Table{1,ifile}(mTable,2)-Table{1,ifile}(1,2);
Frequency{1,ifile}=(mTable/T)*1000;
%Amplitude
meanAmplitude{1,ifile}=mean(Amplitude);
%RiseTime
meanRiseTime{1,ifile}=mean(RiseTime);
%Amplitude CV
AmplitudeCV{1,ifile}=std(Amplitude)/mean(Amplitude);
%IEI
EachTime=Table{1,ifile}(:,2);
IEI{1,ifile}=[0;diff(EachTime)];
%Add Item Name
ItemName{1,ifile}={'number','time(ms)','amplitude','rise(ms)','decay(ms)','area','baseline','noise','group','channel','10-90rise','halfwidth','rise50','peakdir','burst#','burstE#','10-90slope','reltime','IEI'};
CellTable{1,ifile}=[arrayfun(@(x) mat2cell(x),Table{1,ifile}),arrayfun(@(x) mat2cell(x),IEI{1,ifile})];
ZFinalTable{1,ifile}=[ItemName{1,ifile};CellTable{1,ifile}];
%Output Table
FileName{1,ifile}=strrep(ASCname{1,ifile},'.ASC','.xlsx');
xlswrite(OutputName{1,ifile},ZFinalTable{1,ifile});
end
TotalItemName={'Filename','Frequency','Amplitude','RiseTime(ms)','CV_Amplitude','corfficientcorrelation'}';
FinalData=[FileName;Frequency;meanAmplitude;meanRiseTime;AmplitudeCV;corfficientcorrelation];
ZResult=[TotalItemName,FinalData];
xlswrite('All_Result',ZResult);
