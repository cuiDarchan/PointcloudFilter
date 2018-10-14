function Bb = Matrixproperty(C,Bb)
i=length(C);
%% Bb的属性，新添加句
Bb.selection(1:i,:)=1;                              %以下为改写句，把Bb数据的属性值，缩短，原来的B.las索引值1000W
Bb.intensity=ones(i,1,'uint16');                    %换一种数据格式写法
Bb.bits=ones(i,1,'uint8');
Bb.classification=ones(i,1,'uint8');
Bb.user_data=ones(i,1,'uint8');
Bb.scan_angle=ones(i,1,'int8');
Bb.point_source_id=ones(i,1,'uint16');
Bb.gps_time=ones(i,1,'double');
%% 原句
%Bb.selection=Bb.selection(1:i,:);                   %原句
Bb.header.number_of_point_records=i;
Bb.header.number_of_points_by_return(1,1)=i;
Bb.header.max_x=max(C(:,1));
Bb.header.max_y=max(C(:,2));
Bb.header.max_z=max(C(:,3));
Bb.header.min_x=min(C(:,1));
Bb.header.min_y=min(C(:,2));
Bb.header.min_z=min(C(:,3));
Bb.x=C(:,1);
Bb.y=C(:,2);
Bb.z=C(:,3);
%write_las(Bb, 'Bb-denoise4.las', Bb.header.version_major, Bb.header.version_minor, Bb.header.point_data_format)
end



