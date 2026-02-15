This repository contains two image-processing programs:

# A C++ Implementation of BM3D Denoising (Block-Matching 3D Filtering) (BM3D.cpp)
A classic state-of-the-art denoising pipeline that groups similar image patches, applies collaborative filtering in a 3D transform domain, and aggregates overlapping estimates. 
Referred to An Analysis and Implementation of the BM3D Image Denoising Method by Marc Lebrun (2012). Implemented as a two-stage BM3D pipeline:

Step 1: Hard-thresholding collaborative filtering
Step 2: Wiener collaborative filtering

Simply run ./bm3d <file-name> to process an image.

Pre-processed image:

![Alt text](/code/billevans.png?raw=true "BillPre")


Processed image:
![Alt text](/code/finalBill.png?raw=true "BillFinal")

Processed image: 
<img width="1546" height="1068" alt="finalBill" src="https://github.com/user-attachments/assets/cf76fb73-b50c-40da-8c0f-5d488ac97463" />



# A C++ Image Signal Processing Pipeline for RAW Sensor Processing (image.cpp)
A custom image signal processing pipeline designed for processing raw Bayer sensor dumps stored in a specific binary layout (fixed header offset + 16-bit pixel array). It reads raw pixels into a 2D buffer, performs black / white balance, demosaicking (color interpolation), color correction, gamma, tone mapping, edge enhancement, and exports to BMP.

Simply run ./image <file-name> to process an image (Note that images must be of the correct format)

Output image samples:

[DSC00475.SRF.clear.bmp](https://github.com/user-attachments/files/25320109/DSC00475.SRF.clear.bmp)
[DSC00474.SRF.clear.bmp](https://github.com/user-attachments/files/25320108/DSC00474.SRF.clear.bmp)
[DSC00476.SRF.clear.bmp](https://github.com/user-attachments/files/25320110/DSC00476.SRF.clear.bmp)
[DSC00479.SRF.clear.bmp](https://github.com/user-attachments/files/25320111/DSC00479.SRF.clear.bmp)
[DSC00488.SRF.clear.bmp](https://github.com/user-attachments/files/25320112/DSC00488.SRF.clear.bmp)
[DSC00522.SRF.clear.bmp](https://github.com/user-attachments/files/25320113/DSC00522.SRF.clear.bmp)
[DSC00524.SRF.clear.bmp](https://github.com/user-attachments/files/25320115/DSC00524.SRF.clear.bmp)
