Project: 
1.Produce_ground ( main )
2.Groundpoints_classification
3.lasdata
4.Matrixproperty���޸�ͷ�ļ���
5.Reference.las

���ܣ�
�Ի���lidar��las��ʽ�������ݽ����˲���������������ǵ���㡣

ʹ��˵����
1.��Produce_ground�ļ��У��޸ĵ�һ�� B=lasdata('��Ŀ¼����Ҫ�˲���las�ļ�')��
2.��Produce_ground�ļ��У��޸����һ�� 
   write_las(Bb,'��Ŀ¼���˲����las�ļ���',Bb.header.version_major, Bb.header.version_minor, Bb.header.point_data_format);
3.ִ��Produce_ground ���������ɡ�
**˵��ǰ�����У��ļ�����Ҫ���� ".las"��׺

���ԣ�
�����ṩ��һ��test.las�Ĳ������ݣ�������Ϊtest_ground.las���ɽ����������ļ�����envi lidar�в鿴Ч����




