This repository contains two image-processing programs:

# A C++ Implementation of BM3D Denoising (Block-Matching 3D Filtering) (BM3D.cpp)
A classic state-of-the-art denoising pipeline that groups similar image patches, applies collaborative filtering in a 3D transform domain, and aggregates overlapping estimates. Implemented as a two-stage BM3D pipeline:

Step 1: Hard-thresholding collaborative filtering
Step 2: Wiener collaborative filtering

# A C++ Image Signal Processing Pipeline for RAW Sensor Processing (image.cpp)
A custom image signal processing pipeline designed for processing raw Bayer sensor dumps stored in a specific binary layout (fixed header offset + 16-bit pixel array). It reads raw pixels into a 2D buffer, performs black / white balance, demosaicking (color interpolation), color correction, gamma, tone mapping, edge enhancement, and exports to BMP.

Output image samples:
[DSC00544.SRF.clear.bmp](https://github.com/user-a[DSC00524.SRF.clear.bmp](https://github.com/user-attachments/files/25320098/DSC00524.SRF.clear.bmp)
[DSC00522.SRF.clear.bmp](https://github.com/user-attachments/files/25320097/DSC00522.SRF.clear.bmp)
[DSC00488.SRF.clear.bmp](https://github.com/user-attachments/files/25320096/DSC00488.SRF.clear.bmp)
[DSC00479.SRF.clear.bmp](https://github.com/user-attachments/files/25320095/DSC00479.SRF.clear.bmp)
[DSC00476.SRF.clear.bmp](https://github.com/user-attachments/files/25320094/DSC00476.SRF.clear.bmp)
[DSC00475.SRF.clear.bmp](https://github.com/user-attachments/files/25320092/DSC00475.SRF.clear.bmp)
[DSC00474.SRF.clear.bmp](https://github.com/user-attachments/files/25320091/DSC00474.SRF.clear.bmp)
ttachments/files/25320064/DSC00544.SRF.clear.bmp)

