This repository contains two image-processing programs:

# A C++ Implementation of BM3D Denoising (Block-Matching 3D Filtering) (BM3D.cpp)
A classic state-of-the-art denoising pipeline that groups similar image patches, applies collaborative filtering in a 3D transform domain, and aggregates overlapping estimates. Implemented as a two-stage BM3D pipeline:

Step 1: Hard-thresholding collaborative filtering
Step 2: Wiener collaborative filtering

# A C++ Image Signal Processing Pipeline for RAW Sensor Processing (image.cpp)
A custom image signal processing pipeline designed for processing raw Bayer sensor dumps stored in a specific binary layout (fixed header offset + 16-bit pixel array). It reads raw pixels into a 2D buffer, performs black / white balance, demosaicking (color interpolation), color correction, gamma, tone mapping, edge enhancement, and exports to BMP.

Output image samples:

[DSC00475.SRF.clear.bmp](https://github.com/user-attachments/files/25320109/DSC00475.SRF.clear.bmp)
[DSC00474.SRF.clear.bmp](https://github.com/user-attachments/files/25320108/DSC00474.SRF.clear.bmp)
[DSC00476.SRF.clear.bmp](https://github.com/user-attachments/files/25320110/DSC00476.SRF.clear.bmp)
[DSC00479.SRF.clear.bmp](https://github.com/user-attachments/files/25320111/DSC00479.SRF.clear.bmp)
[DSC00488.SRF.clear.bmp](https://github.com/user-attachments/files/25320112/DSC00488.SRF.clear.bmp)
[DSC00522.SRF.clear.bmp](https://github.com/user-attachments/files/25320113/DSC00522.SRF.clear.bmp)
[DSC00524.SRF.clear.bmp](https://github.com/user-attachments/files/25320115/DSC00524.SRF.clear.bmp)
