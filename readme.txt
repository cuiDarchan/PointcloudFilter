Project: 
1.Produce_ground ( main )
2.Groundpoints_classification
3.lasdata
4.Matrixproperty（修改头文件）
5.Reference.las

功能：
对机载lidar的las格式点云数据进行滤波，即分离地面点与非地面点。

使用说明：
1.在Produce_ground文件中，修改第一行 B=lasdata('本目录中需要滤波的las文件')；
2.在Produce_ground文件中，修改最后一行 
   write_las(Bb,'本目录中滤波后的las文件名',Bb.header.version_major, Bb.header.version_minor, Bb.header.point_data_format);
3.执行Produce_ground 主函数即可。
**说明前两条中，文件名需要加上 ".las"后缀

测试：
本例提供了一个test.las的测试数据，处理结果为test_ground.las，可将上述两个文件放入envi lidar中查看效果。




