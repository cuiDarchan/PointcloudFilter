function [groundpoints,nogroundpoints]=Groundpoints_classification(Ave_Num,Ave_m,B,threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%groundpoints  ��ȡ�ĵ������
%nogroundpoints�ǵ������
%Ave_Num       ÿ��������ϣ���ĵ��Ƹ���
%ave_m:        �����ܶ�
%B:            BΪ����ĵ��ƾ���txt��һ��Ϊn*4��ʽ
%threshold     ��ֵ��������ȡ2
%Rmin          �ٶ��̱߳仯����ķ�Χ����Сֵ
%Rmax          �ٶ��̱߳仯����ķ�Χ�����ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��������
%Ave_Num=20;
ave=fix(sqrt((Ave_Num/Ave_m)))+1;
% chajz=zeros(length(B),1);
% tt=1;
d=1;                                           %���������ֵ
nd=1;                                          %�ǵ��������ֵ
h_threshlod=3;                                 %�ٶ�������ڸ߳�z��ֵ����ֵ�����������ֱ���ж�Ϊ�ǵ����
delt_threshold=3;                              %�ٶ�������е�С�������ᳬ���������������ֵ3m����
zg=0;                                          %z�������Ĺ�ֵ
%% ��ȡ16����������͵㣬�������Qm�����棩
% ��������B����Ϊ����n*5��ʽ,�ֱ�Ϊx,y,z,�Լ���ά��������ֵ
len=length(B);
I=zeros(len,6);
for i=1:len
    I(i,1)=B(i,1);
    I(i,2)=B(i,2);
    I(i,3)=B(i,3);
    %I(i,6)=B(i,4);
end
% д���Ԫ����
xmin=min(I(:,1));
xmax=max(I(:,1));
ymin=min(I(:,2));
ymax=max(I(:,2));
ZMIN=min(I(:,3));                                %ȫ�����ֵ
for t=1:len
    I(t,4)=fix((I(t,1)-xmin)/ave)+1;
    I(t,5)=fix((I(t,2)-ymin)/ave)+1;
end
%% �����������ڲ����װ����Ӧ����Ԫ���ڣ���������
M=max(I(:,4));
N=max(I(:,5));
netcell=cell(M,N);                               %��M,��N����ά�����߶�ȷ��
count=zeros(M,N);
for i=1:len
    [p,~]=size(netcell{I(i,4),I(i,5)});          %��ʼʱ��pΪ0�������� ����д���Ԫ��
    netcell{I(i,4),I(i,5)}(p+1,1)=I(i,1);        %langcell{}��ָ���ð�Ԫ��Ԫ��,��()������
    netcell{I(i,4),I(i,5)}(p+1,2)=I(i,2);
    netcell{I(i,4),I(i,5)}(p+1,3)=I(i,3);
    count(I(i,4),I(i,5))=count(I(i,4),I(i,5))+1;
end
%%  ��ѡ6����ϵ�ı���
%�������
Nh=6;                                            %NhΪ��ϵ����
groundpoints=zeros(len,5);                       %Ϊ���������ڴ�ռ�
nogroundpoints=zeros(len,5);
theshold_Num=4;                                  %��Ѱ�����м��4�����ֵ
Qm=zeros(Nh,3);                                  %6��������Z��һ�͵�ֵ�������������
Rmin=6*ave;
Rmax=10*ave;
Rf=(fix(Rmin/ave)+1)*ave;
Rl=(fix(Rmax/ave))*ave;
R_num=(Rl-Rf)/ave+1;
RMatrix=zeros(R_num,4);
R_i=1;
for R=Rf:+ave:Rl
    RMatrix(R_i,1)=R;
    RMatrix(R_i,2)=mod(M*ave,R);
    RMatrix(R_i,3)=mod(N*ave,R);
    RMatrix(R_i,4)=RMatrix(R_i,2)+RMatrix(R_i,3);
    R_i=R_i+1;
end
R_min=min(RMatrix(:,4));
[row,~]=find(RMatrix(:,4)==R_min);
avec=RMatrix(row(1),1);
fg=avec/ave;                                  %������зָ��С��������
MM=fix(M*ave/avec);                           %������߶��зָ���
NN=fix(N*ave/avec);
RefNum=cell(MM+1,NN+1);                       %���c0-c5����������Ԫ��
Zmi=zeros(MM+1,NN+1);                         %Zmi��������Z�߳��ĵ㣬�Ż��ٶ�
Pos=zeros(MM+1,NN+1);                         %ZpcΪ��Ӧ���������ڵ���������
var_th=5;                                     %���ݲ�����ֵ
hz_th=1;                                      %ѡȡ��������ֵ
yc_Matrix=ones(MM,NN);                        %�쳣������󣬼�¼��������
yc_Matrix2=ones(MM,NN);                       %�쳣�������2����¼��������
X=zeros(Nh,Nh);                               %����������
JG=zeros(Nh,1);                               %�������
%% �����������Ѱ����
for ii=1:MM
    for jj=1:NN
        while 1
            Num=fix(fg*rand(1,1))+1;            
            if Num>0.8*fg
                Num=Num-4;
            elseif Num<0.2*fg
                Num=Num+4;
            end
            %�쳣���������Ԫ��Ϊ�յ����
            if isempty(netcell{Num+(ii-1)*fg,Num+(jj-1)*fg})==1||isempty(netcell{Num+(ii-1)*fg,Num+1+(jj-1)*fg})==1||isempty(netcell{Num+(ii-1)*fg,Num+2+(jj-1)*fg})==1 ...
                    ||isempty(netcell{Num+2+(ii-1)*fg,Num+(jj-1)*fg})==1||isempty(netcell{Num+2+(ii-1)*fg,Num+1+(jj-1)*fg})==1||isempty(netcell{Num+2+(ii-1)*fg,Num+2+(jj-1)*fg})==1
                yc_Matrix(ii,jj)=yc_Matrix(ii,jj)+1;
                if yc_Matrix(ii,jj)>=10        %�������10�λ�δ�ҵ�Num������һ�����е�Qm
                    break;
                end
                continue;
            end
            Pos(ii,jj)=Num;                    %��Pos�м�¼�����6������λ��Num
            %�ҵ�6����Сֵ��
            x=1;
            for u=Num:+2:Num+2
                for l=Num:Num+2
                    netcell{u+(ii-1)*fg,l+(jj-1)*fg}=sortrows(netcell{u+(ii-1)*fg,l+(jj-1)*fg},3);
                    Qm(x,:)=netcell{u+(ii-1)*fg,l+(jj-1)*fg}(1,:);
                    x=x+1;
                end
            end
            %�б���
            if jj==1&&ii==1
                bl1=mean(Qm(:,3))-(ZMIN+10);
            elseif jj~=1
                bl1=mean(Qm(:,3))-(Zmi(ii,jj-1)+3);
            elseif jj==1&&ii~=1
                bl1=mean(Qm(:,3))-(Zmi(ii-1,jj)+3);
            end
            bl2=var(Qm(:,3))-var_th;
            yc_Matrix2(ii,jj)=yc_Matrix2(ii,jj)+1;        %�쳣�������
            if yc_Matrix2(ii,jj)>=10                      %����10��ֹͣ
                break;
            end
            if (bl1<0&&bl2<0)
                break;
            end
        end
        %��ʼ�����6����
        Center=mean(Qm);                                  %���Qm���������
        xp0=Center(1);yp0=Center(2);zp0=Center(3);
        for q=1:Nh
            X(q,1)=1;
            X(q,2)=Qm(q,1)-xp0;
            X(q,3)=Qm(q,2)-yp0;
            X(q,4)=(Qm(q,1)-xp0).^2;
            X(q,5)=(Qm(q,1)-xp0)*(Qm(q,2)-yp0);
            X(q,6)=(Qm(q,2)-yp0).^2;
            JG(q,1)=Qm(q,3)-zp0;
        end
        A=X\JG;              
        RefNum{ii,jj}=A;                                  %����õĲ�����ŵ���Ӧ������
        Zmi(ii,jj)=zp0;                                   %ZmiΪ����߳����ģ���������Z�߳��ĵ㣬�Ż��ٶ�
    end
end
%3����������ֵRefNum��Zmi��Pos
RefNum(ii+1,:)=RefNum(ii,:);RefNum(:,jj+1)=RefNum(:,jj);
Zmi(ii+1,:)=Zmi(ii,:);Zmi(:,jj+1)=Zmi(:,jj);
Pos(ii+1,jj)=1; Pos(ii+1,jj+1)=1; Pos(ii,jj+1)=1;         %���߽紦��Pos��Num���ó�1
%% ��������б�
%��һ��������С�����ڵĶ�ά�������꣨x��y��,������Position�У����������ִ�����������ֵ
Position=cell(M,N);
for i=1:M
    for j=1:N
        di=fix((i-1)/fg)+1;                               %������µ���������di��dj
        dj=fix((j-1)/fg)+1;
        Position{i,j}(1,1)=xmin+(i-0.5)*ave;
        Position{i,j}(1,2)=ymin+(j-0.5)*ave;
        %�쳣����С�����ڵ�������6ʱ��,�����������Ϊ0
        [g,~]=size(netcell{i,j});
        if g<Nh
            Position{i,j}(1,3)=0;
            continue;
        end
        netcell{i,j}=sortrows(netcell{i,j},3);
        Position{i,j}(1,3)=mean(netcell{i,j}(3:6,3));
        if Position{i,j}(1,3)>Zmi(di,dj)+delt_threshold   %�쳣�㴦�������Ƿ��ݻ��߶���
            Position{i,j}(1,3)=0;
        end
    end
end
%�ڶ�����������ΧԪ��������쳣Ԫ���ĸ߳�ֵ���滻0�����ܿ�߽�δ������
for i=2:M-1
    for j=2:N-1
        if  Position{i,j}(1,3)==0
        Ng=9;
        sum=0;
        for k=-1:1
            for l=-1:1
               if Position{i+k,j+l}(1,3)==0
                   Ng=Ng-1;
               else
                   sum=sum+Position{i+k,j+l}(1,3);
               end                          
            end
        end
        %�쳣����
        if Ng==0
            continue;
        else
        Position{i,j}(1,3)=sum/Ng;
        end
        end
    end
end
%�߽����,�ٶ�4�ǲ���Ҫ����
 %�����һ��
for j=2:N-1
    if Position{1,j}(1,3)==0
        sum=0;
        Ng=6;
        for k=0:1
            for l=-1:1
                if Position{1+k,j+l}(1,3)==0
                    Ng=Ng-1;
                else
                    sum=sum+Position{1+k,j+l}(1,3);
                end
            end
        end
        %�쳣����,�����ò���
        if Ng==0                        
            continue;
        else
        Position{1,j}(1,3)=sum/Ng;
        end
    end
end
 %��������һ��
for j=2:N-1
    if Position{M,j}(1,3)==0
        sum=0;
        Ng=6;
        for k=-1:0
            for l=-1:1
                if Position{M+k,j+l}(1,3)==0
                    Ng=Ng-1;
                else
                    sum=sum+Position{M+k,j+l}(1,3);
                end
            end
        end
        %�쳣����,�����ò���
        if Ng==0                        
            continue;
        else
        Position{M,j}(1,3)=sum/Ng;
        end
    end
end
 %���������
for i=2:N-1
    if Position{i,1}(1,3)==0
        sum=0;
        Ng=6;
        for k=-1:1
            for l=0:1
                if Position{i+k,1+l}(1,3)==0
                    Ng=Ng-1;
                else
                    sum=sum+Position{i+k,1+l}(1,3);
                end
            end
        end
        %�쳣����,�����ò���
        if Ng==0                        
            continue;
        else
        Position{i,1}(1,3)=sum/Ng;
        end
    end
end
 %�������Ҳ�
for i=2:N-1
    if Position{i,N}(1,3)==0
        sum=0;
        Ng=6;
        for k=-1:1
            for l=-1:0
                if Position{i+k,N+l}(1,3)==0
                    Ng=Ng-1;
                else
                    sum=sum+Position{i+k,N+l}(1,3);
                end
            end
        end
        %�쳣����,�����ò���
        if Ng==0                        
            continue;
        else
        Position{i,N}(1,3)=sum/Ng;
        end
    end
end
%���飬�Ƿ����0ֵ
% for i=1:M
%     for j=1:N
%         if Position{i,j}(1,3)==0
%            AA=i;
%            AA=j;
%         end
%     end
% end

%����������һ���newpoint
newpoint=zeros(1,3);
for i=1:M
    for j=1:N
        ii1=fix((i-1)/fg)+1;
        jj1=fix((j-1)/fg)+1;
        [p,~]=size(netcell{i,j});
        for k=1:p
            newpoint=netcell{i,j}(k,:);
            if newpoint(3)<Zmi(ii1,jj1)+h_threshlod;
                %�������
                xp=newpoint(1);yp=newpoint(2);zp=newpoint(3);
                zg=A(2)*(xp-Position{i,j}(1))+A(3)*(yp-Position{i,j}(2))+A(4)*(xp-Position{i,j}(1)).^2+A(5)*(xp-Position{i,j}(1)).*(yp-Position{i,j}(2))+A(6)*(yp-Position{i,j}(2)).^2+Position{i,j}(3);
                cha=double(zg-zp);                      %chaΪ�ӽ�0����
                %�鿴cha�����ֵ��Ѱ�Ҿ�����ֵ
%                 chajz(tt,1)=cha;
%                 tt=tt+1;
                if abs(cha)<threshold
                    groundpoints(d,1:3)=newpoint;
                    d=d+1;
                else
                    nogroundpoints(nd,1:3)=newpoint;
                    nd=nd+1;
                end
            else
                nogroundpoints(nd,1:3)=newpoint;
                nd=nd+1;
            end
            
        end
        
    end
end
groundpoints=groundpoints(1:d-1,:);
nogroundpoints=nogroundpoints(1:nd-1,:);
end
