def dcm2jpg(path, out_path):
    import SimpleITK as sitk
    import numpy as np
    import cv2
    import os


    import time
    from PIL import Image
    path = "E:\\Mask_RCNN_master\\UPENN-GBM-DWI\\"
    out_path = "E:\\Mask_RCNN_master\\UPENN-GBM-DWI-JPG\\"
    filename = os.listdir(path)
    for file in filename:
        path1 = path + file + '/'
        filename_1 = os.listdir(path1)
        for file1 in filename_1:
            outputpath = out_path + file + "/" + file1 + "/"
            os.makedirs(outputpath)
            shibai_path = path1 + file1 + "/"
            c = []
            names = os.listdir(shibai_path)  # 路径
    # 将文件夹中的文件名称与后边的 .dcm分开

            for name in names:
                index = name.rfind('.')
                name = name[:index]
                c.append(name)

            for i in c:
                document = shibai_path + i + ".dcm"
                countname = i
                countfullname = countname + '.jpg'
                output_jpg_path = os.path.join(outputpath, countfullname)


                def convert_from_dicom_to_jpg(img, low_window, high_window, save_path):
                    lungwin = np.array([low_window * 1., high_window * 1.])
                    newimg = (img - lungwin[0]) / (lungwin[1] - lungwin[0])
                    newimg = (newimg * 255).astype('uint8')
                    newimg = cv2.resize(newimg, (256, 256), interpolation=cv2.INTER_CUBIC)
                    cv2.imwrite(save_path, newimg, [int(cv2.IMWRITE_JPEG_QUALITY), 100])


                if __name__ == '__main__':
                    ds_array = sitk.ReadImage(document)
                    img_array = sitk.GetArrayFromImage(ds_array)
                    shape = img_array.shape  # name.shape
                    img_array = np.reshape(img_array, (shape[1], shape[2]))
                    high = np.max(img_array)
                    low = np.min(img_array)
                    convert_from_dicom_to_jpg(img_array, low, high, output_jpg_path)
                    print('FINISHED')

dcm2jpg("E:\\Mask_RCNN_master\\UPENN-GBM-DWI\\", "E:\\Mask_RCNN_master\\UPENN-GBM-DWI-JPG\\")