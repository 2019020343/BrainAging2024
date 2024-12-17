import os
from shutil import copyfile

for root, files, names in os.walk(r'E:\pancancer\Awesome-Backbones-main\pancancer_photos_gbm\result\TUMOR'):  # 改成你自己的json文件夹所在的目录
    for name in names:
        sample_name = name.split("_")[0]
        file_path = r'E:\Mask_RCNN_master\TCGA-GBM-T1-JPG'  # 目标路径
        file_dir = os.path.join(file_path, sample_name)
        new_file_name = os.path.join(file_dir, name)

        tar_root = r'E:\Mask_RCNN_master\TCGAGBM_Resnet_all_T_original'  # 目标路径
            #tar_root = r'E:\Mask_RCNN_master\HGG_train_TUMOR\label1'  # 目标路径
        tar_dir = os.path.join(tar_root, sample_name)
        if not os.path.exists(tar_dir):
            os.makedirs(tar_dir)
        tar_file = os.path.join(tar_dir, name)
        if not os.path.exists(new_file_name):
            print(f"文件 {new_file_name} 不存在，跳过。")
        else:
            copyfile(new_file_name, tar_file)