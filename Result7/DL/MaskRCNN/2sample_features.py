# -*- coding: utf-8 -*-


import os  # 读取txt文件所需要的包
import linecache  # 读取指定行函数linecache.getline(file_ob, line_num)所在的包

# 1. 将要处理的文件数据txt放入pyradiomics文件下：
root_path = 'E:\\Mask_RCNN_master\\pyradiomics' # 读取的批量txt所在的总文件夹路径
features_root = os.path.join(root_path, "T1_features") #读取的批量txt所在的特征文件夹路径
# features_root = os.path.join(root_path, "test") #读取的批量txt所在的特征文件夹路径
masked_root = os.path.join(root_path, "T1_masked") #读取的批量txt所在的掩膜文件夹路径
sample_names = os.listdir(masked_root)  # 读取pyradiomics文件夹下的样本文件名
for sample in sample_names:
    sample_path = os.path.join(masked_root, sample)
    file_names = os.listdir(sample_path)
    if file_names.__len__()>0:
      # 读取pyradiomics文件夹下的txt文件名
        file_ob_list = []  # 定义一个列表，用来存放刚才读取的txt文件名
        for file_name in file_names:
            fileob = features_root + '\\' + sample + '\\' + file_name.split(".")[0] + ".txt"  ##文件夹路径加上\\ 再加上具体要读的的txt的文件名就定位到了这个txt
            file_ob_list.append(fileob)  # 将路径追加到列表中存储
            # 这里添加路径，方便后续linecache.getline()进行数据提取
        print(file_ob_list)  # 打印这个列表的内容到显示屏，不想显示的话可以去掉这句

        ldata = []  # 收集所有行数据
        data = []  # 收集每一行数据，每次循环后，都要清空

        # file_ob_list是所有文件（比如10个txt）对象组成的列表，for用来循环读取每一个文件，读取一个文件的方式是一行行读入，
        # 每次循环一次for就读取所有文件的某一行，因为这一行的第一列都是特征名称，都是一样的

        line_num = 23  # 从txt的第23行开始读入
        total_line = len(open(file_ob_list[0]).readlines())  # 计算一个txt中有多少行
        while line_num <= total_line:  # 只有读完的行数小于等于总行数时才再读下一行，否则结束读取
            for file_ob in file_ob_list:  # 按顺序循环读取所有文件
                line = linecache.getline(file_ob, line_num)  # 读取这个文件的第line_num行
                line = line.strip()  # 去掉这一行最后一个字符\n，即换行符
                if line is None or len(line) == 0:
                    break
                fields = line.split('\t')  # 将这一行划分为两列，存放到列表中，fields是这样的['original_firstorder_10Percentile','63.60000000000002']
                prob = fields[1]  # fields[0]是'original_firstorder_10Percentile'   fields[1]是'63.60000000000002'

                # 这个if部分只是将特征量的10位小数点压缩到4位，其实可以去掉这个处理
                if fields[1] != "VALUE":  # 基因表达量不是NA也就是为数字时，才对它进行小数点的减少处理
                    prob = float(fields[1])  # 将字符形式的数字如'0.0'强制转化为浮点型（带小数点的数字）数字0.0
                    prob = '%.4f' % prob  # 只保留小数点后面的4位小数

                if file_ob == file_ob_list[0]:  # 如果读的是第一个txt文件，则将读进去的第一列基因名和第二列表达量
                    data = [fields[0], prob]  # 都加入到列表中 data= ['original_firstorder_10Percentile','63.6000']
                else:  # 如果读进去的不是第一个文件，则跳到else执行，第一列不要，
                    data.append(prob)  # 只将第二列表达量追加到之前的二维数组后面,假如这时读的是第二个文件的第一行
                    # 此时fields为['original_firstorder_10Percentile','19.9000']，则data=['original_firstorder_10Percentile','63.6000'，'19.9000']
            line_num = line_num + 1  # 行数加1，好接着读取每一个文件的第二行 (每个文件逐行读入，并存入)
            ldata.append(data)  # 将存放了所有txt的第一行数据的data，放到一个新的列表中保存，这时ldata就是一个二维列表，ldata=[['original_firstorder_10Percentile','63.6000'，'19.9000'],[...],...]
            # 用来存放所有的（108行（前22行是软件库版本等注释信息， 后86行是特征，），就是所有的特征名）行数
            data = []  # 清空data用来存放所有txt的下一行
        out_path = root_path + "\\" + "T1_sample_feature" + "\\"  + sample + "_result.txt"
        f = open(out_path, "w+")  # 创建存放数据的文件，目前是空的，需要进一步写入

    #    ldata.pop(0)  # 删除数据自带的header文件

        # 将数据加如header
        file_names_new1 = []
        for col_name in file_names:
            col_name = col_name.strip("_sample_table.txt")
            file_names_new1.append(col_name)
        file_names_new2 = ["Feature", ] + file_names_new1

        temp = [file_names_new2, ] + ldata
        # temp.append(ldata) #合并header文件与ldata的value文件

        # 写入数据
        for i, p in enumerate(temp):  # 将数据写入文件,i是enumerate()函数返回的ldata的某个元素p(就是一行数据，如['ENSG242268.2','0.0'，'0.10']从第一个开始)开始的序号（0,1,2等）
            for j, q in enumerate(p):  # 读取p（如['ENSG242268.2','0.0'，'0.10']）中的每一个元素
                if j == len(temp[i]) - 1:
                    f.write(q)  # 将这个元素写到txt中，每写一个加入一个“\t”（它代表excel中的一根竖线）
                else:
                    f.write(q + "\t")  # 将这个元素写到txt中，每写一个加入一个“\t”（它代表excel中的一根竖线）
            print(i)  # 显示一下打印到了第多少行
            f.write("\n")  # 每写完一行，就写入一个换行符“\n”， 好使接下来的数据写入到第二行

        f.close()  # 操作完一个文件后应该将其关闭