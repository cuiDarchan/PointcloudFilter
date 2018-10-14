function [groundpoints,nogroundpoints]=Groundpoints_classification(Ave_Num,Ave_m,B,threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%groundpoints  获取的地面点云
%nogroundpoints非地面点云
%Ave_Num       每个格网内希望的点云个数
%ave_m:        点云密度
%B:            B为输入的点云矩阵，txt下一般为n*4格式
%threshold     阈值，经验上取2
%Rmin          假定高程变化不大的范围，最小值
%Rmax          假定高程变化不大的范围，最大值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 变量定义
%Ave_Num=20;
ave=fix(sqrt((Ave_Num/Ave_m)))+1;
% chajz=zeros(length(B),1);
% tt=1;
d=1;                                           %地面点索引值
nd=1;                                          %非地面点索引值
h_threshlod=3;                                 %假定大格网内高出z均值的阈值，如果超出，直接判断为非地面点
delt_threshold=3;                              %假定大格网中的小格网不会超过大格网参数中心值3m以上
zg=0;                                          %z相对坐标的估值
%% 获取16个格网内最低点，存入矩阵Qm（曲面）
% 加载数据B，变为正常n*5格式,分别为x,y,z,以及二维格网坐标值
len=length(B);
I=zeros(len,6);
for i=1:len
    I(i,1)=B(i,1);
    I(i,2)=B(i,2);
    I(i,3)=B(i,3);
    %I(i,6)=B(i,4);
end
% 写入胞元数据
xmin=min(I(:,1));
xmax=max(I(:,1));
ymin=min(I(:,2));
ymax=max(I(:,2));
ZMIN=min(I(:,3));                                %全域最低值
for t=1:len
    I(t,4)=fix((I(t,1)-xmin)/ave)+1;
    I(t,5)=fix((I(t,2)-ymin)/ave)+1;
end
%% 将各个格网内部点封装到相应格网元胞内，方便索引
M=max(I(:,4));
N=max(I(:,5));
netcell=cell(M,N);                               %长M,宽N，二维格网尺度确定
count=zeros(M,N);
for i=1:len
    [p,~]=size(netcell{I(i,4),I(i,5)});          %初始时，p为0，将长宽 数据写入胞元内
    netcell{I(i,4),I(i,5)}(p+1,1)=I(i,1);        %langcell{}是指引用胞元内元素,与()有区别
    netcell{I(i,4),I(i,5)}(p+1,2)=I(i,2);
    netcell{I(i,4),I(i,5)}(p+1,3)=I(i,3);
    count(I(i,4),I(i,5))=count(I(i,4),I(i,5))+1;
end
%%  挑选6个拟合点的变量
%定义变量
Nh=6;                                            %Nh为拟合点个数
groundpoints=zeros(len,5);                       %为地面点分配内存空间
nogroundpoints=zeros(len,5);
theshold_Num=4;                                  %找寻曲面中间的4个最低值
Qm=zeros(Nh,3);                                  %6个格网内Z第一低的值，用于拟合曲面
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
fg=avec/ave;                                  %大格网中分割的小格网个数
MM=fix(M*ave/avec);                           %格网大尺度切分个数
NN=fix(N*ave/avec);
RefNum=cell(MM+1,NN+1);                       %存放c0-c5六个参数的元胞
Zmi=zeros(MM+1,NN+1);                         %Zmi用来过滤Z高出的点，优化速度
Pos=zeros(MM+1,NN+1);                         %Zpc为对应各个网格内的曲面中心
var_th=5;                                     %数据波动阈值
hz_th=1;                                      %选取随机点的阈值
yc_Matrix=ones(MM,NN);                        %异常处理矩阵，记录迭代次数
yc_Matrix2=ones(MM,NN);                       %异常处理矩阵2，记录迭代次数
X=zeros(Nh,Nh);                               %求解参数矩阵
JG=zeros(Nh,1);                               %结果矩阵
%% 遍历大格网找寻参数
for ii=1:MM
    for jj=1:NN
        while 1
            Num=fix(fg*rand(1,1))+1;            
            if Num>0.8*fg
                Num=Num-4;
            elseif Num<0.2*fg
                Num=Num+4;
            end
            %异常处理，如果胞元内为空的情况
            if isempty(netcell{Num+(ii-1)*fg,Num+(jj-1)*fg})==1||isempty(netcell{Num+(ii-1)*fg,Num+1+(jj-1)*fg})==1||isempty(netcell{Num+(ii-1)*fg,Num+2+(jj-1)*fg})==1 ...
                    ||isempty(netcell{Num+2+(ii-1)*fg,Num+(jj-1)*fg})==1||isempty(netcell{Num+2+(ii-1)*fg,Num+1+(jj-1)*fg})==1||isempty(netcell{Num+2+(ii-1)*fg,Num+2+(jj-1)*fg})==1
                yc_Matrix(ii,jj)=yc_Matrix(ii,jj)+1;
                if yc_Matrix(ii,jj)>=10        %如果迭代10次还未找到Num，用上一方格中的Qm
                    break;
                end
                continue;
            end
            Pos(ii,jj)=Num;                    %在Pos中记录计算出6参数的位置Num
            %找到6个最小值点
            x=1;
            for u=Num:+2:Num+2
                for l=Num:Num+2
                    netcell{u+(ii-1)*fg,l+(jj-1)*fg}=sortrows(netcell{u+(ii-1)*fg,l+(jj-1)*fg},3);
                    Qm(x,:)=netcell{u+(ii-1)*fg,l+(jj-1)*fg}(1,:);
                    x=x+1;
                end
            end
            %判别处理
            if jj==1&&ii==1
                bl1=mean(Qm(:,3))-(ZMIN+10);
            elseif jj~=1
                bl1=mean(Qm(:,3))-(Zmi(ii,jj-1)+3);
            elseif jj==1&&ii~=1
                bl1=mean(Qm(:,3))-(Zmi(ii-1,jj)+3);
            end
            bl2=var(Qm(:,3))-var_th;
            yc_Matrix2(ii,jj)=yc_Matrix2(ii,jj)+1;        %异常处理矩阵
            if yc_Matrix2(ii,jj)>=10                      %迭代10次停止
                break;
            end
            if (bl1<0&&bl2<0)
                break;
            end
        end
        %初始化求解6参数
        Center=mean(Qm);                                  %求解Qm矩阵的中心
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
        RefNum{ii,jj}=A;                                  %将求得的参数存放到对应坐标下
        Zmi(ii,jj)=zp0;                                   %Zmi为曲面高程中心，用来过滤Z高出的点，优化速度
    end
end
%3个变量矩阵赋值RefNum，Zmi，Pos
RefNum(ii+1,:)=RefNum(ii,:);RefNum(:,jj+1)=RefNum(:,jj);
Zmi(ii+1,:)=Zmi(ii,:);Zmi(:,jj+1)=Zmi(:,jj);
Pos(ii+1,jj)=1; Pos(ii+1,jj+1)=1; Pos(ii,jj+1)=1;         %将边界处的Pos的Num都置成1
%% 曲面拟合判别
%第一步，计算小格网内的二维中心坐标（x，y）,并存入Position中，第三个数字存入格网内最低值
Position=cell(M,N);
for i=1:M
    for j=1:N
        di=fix((i-1)/fg)+1;                               %大格网下的坐标索引di，dj
        dj=fix((j-1)/fg)+1;
        Position{i,j}(1,1)=xmin+(i-0.5)*ave;
        Position{i,j}(1,2)=ymin+(j-0.5)*ave;
        %异常处理，小格网内点数少于6时候,令第三个数字为0
        [g,~]=size(netcell{i,j});
        if g<Nh
            Position{i,j}(1,3)=0;
            continue;
        end
        netcell{i,j}=sortrows(netcell{i,j},3);
        Position{i,j}(1,3)=mean(netcell{i,j}(3:6,3));
        if Position{i,j}(1,3)>Zmi(di,dj)+delt_threshold   %异常点处理，可能是房屋或者陡坡
            Position{i,j}(1,3)=0;
        end
    end
end
%第二步，利用周围元胞，求出异常元胞的高程值，替换0，四周框边界未做处理
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
        %异常处理
        if Ng==0
            continue;
        else
        Position{i,j}(1,3)=sum/Ng;
        end
        end
    end
end
%边界框处理,假定4角不需要处理
 %处理第一行
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
        %异常处理,估计用不到
        if Ng==0                        
            continue;
        else
        Position{1,j}(1,3)=sum/Ng;
        end
    end
end
 %处理最下一行
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
        %异常处理,估计用不到
        if Ng==0                        
            continue;
        else
        Position{M,j}(1,3)=sum/Ng;
        end
    end
end
 %处理最左侧
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
        %异常处理,估计用不到
        if Ng==0                        
            continue;
        else
        Position{i,1}(1,3)=sum/Ng;
        end
    end
end
 %处理最右侧
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
        %异常处理,估计用不到
        if Ng==0                        
            continue;
        else
        Position{i,N}(1,3)=sum/Ng;
        end
    end
end
%检验，是否存在0值
% for i=1:M
%     for j=1:N
%         if Position{i,j}(1,3)==0
%            AA=i;
%            AA=j;
%         end
%     end
% end

%第三步，逐一添加newpoint
newpoint=zeros(1,3);
for i=1:M
    for j=1:N
        ii1=fix((i-1)/fg)+1;
        jj1=fix((j-1)/fg)+1;
        [p,~]=size(netcell{i,j});
        for k=1:p
            newpoint=netcell{i,j}(k,:);
            if newpoint(3)<Zmi(ii1,jj1)+h_threshlod;
                %计算参数
                xp=newpoint(1);yp=newpoint(2);zp=newpoint(3);
                zg=A(2)*(xp-Position{i,j}(1))+A(3)*(yp-Position{i,j}(2))+A(4)*(xp-Position{i,j}(1)).^2+A(5)*(xp-Position{i,j}(1)).*(yp-Position{i,j}(2))+A(6)*(yp-Position{i,j}(2)).^2+Position{i,j}(3);
                cha=double(zg-zp);                      %cha为接近0的数
                %查看cha矩阵的值，寻找经验阈值
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
