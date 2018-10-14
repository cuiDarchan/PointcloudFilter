B=lasdata('test.las');
[i,~]=size(B.selection);
C=zeros(i,3);
C(1:i,1)=B.x;
C(1:i,2)=B.y;
C(1:i,3)=B.z;
%点云密度
Xmin=min(C(:,1));
Xmax=max(C(:,1));
Ymin=min(C(:,2));
Ymax=max(C(:,2));
pm=i/((Xmax-Xmin).*(Ymax-Ymin));
[groundpoints,nogroundpoints]=Groundpoints_classification(20,pm,C,0.5);
Bb=lasdata('Reference.las');     
Bb=Matrixproperty(groundpoints,Bb);                  %改写头文件
write_las(Bb,'test2.las',Bb.header.version_major, Bb.header.version_minor, Bb.header.point_data_format);