# Deep learning Segmentation Reproducibility Study
**The Reproducibility of Deep Learning-Based Segmentation of the Prostate Gland and Zones on T2-Weighted MR Images**

This is the full script that used to conduct the study.

The study was done at the MR Cancer group at the Norwegian University of Science and Technology (NTNU) in Trondheim, Norway.
https://www.ntnu.edu/isb/mr-cancer

For detailed information about this method, please read our paper: https://www.mdpi.com/2075-4418/11/9/1690

# Note
The provided script was used for research use only.

# How to cite this work
In case of using or refering to this script or study, please cite it as:

Sunoqrot, M.R.S.; Selnæs, K.M.; Sandsmark, E.; Langørgen, S.; Bertilsson, H.; Bathen, T.F.; Elschot, M. The Reproducibility of Deep Learning-Based Segmentation of the Prostate Gland and Zones on T2-Weighted MR Images. Diagnostics 2021, 11, 1690.
https://doi.org/10.3390/diagnostics11091690

# How to use the script
This is a MATLAB® script, the script was written and tested using MATLAB® R2019b.

The file "Master.m" is the main script that contains all the sub functions of the analysis. It also allows to control which fuctions to run or not.

Make sure that all of these files are in the same folder.

**Input:**
You will need to change the paths in the script, mainly the Master file and make sure you prepared the data according to the first section on the analysis script in Master.

**Output:**
 Statistical analysis report with some figures, tables and examples.
  
# Dependency 
This script depend on the followings, which you should make sure that you have correctly installed them on your computer:
1. AutoRef normalization method
  by MR Cancer Group at NTNU Trondheim, Norway https://github.com/ntnu-mr-cancer/AutoRef
  - It is not included in the Dependency foldr.
  - Make sure to install it and normalize the images using it, since you will need them to start the invistigation.
  You can find the method and all the informatio you need at:
  https://github.com/ntnu-mr-cancer/PSQC/tree/master/Dependency/AutoRef
2. JSONLab toolbox (version 1.5):
  by:  Qianqian Fang, Northeastern University, MA, USA.
  - It is included in the Dependency folder.
3.  Python environment with Pyradiomics (V 3.0) and python (3.7).
  Pyradiomics is by: Computational Imaging & Bioinformatics Lab. Harvard Medical School, MA , USA.
  https://pyradiomics.readthedocs.io/en/2.2.0/
4.  Convert3D tool from ITK 
  by ITK-SNAP http://www.itksnap.org
  - You MUST install it to your computer and complie it with your system as descriped: 
    + Install the "CONVERT3D NIGHTLY BUILD" folder from: http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.C3D
    + Follow the guide in the documents to install and compile: http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Convert3D  
5.  elastix toolbox (4.3<=version<=4.7):
  by: Image Sciences Institute, University Medical Center Utrecht, The Netherlands.
  - It is included in the Dependency folder, so no need to download it unless you faced a problem with it, in that case:
    + Download and compile as descriped at: http://elastix.isi.uu.nl/
    + Read and follow the section 1.2 at: http://elastix.isi.uu.nl/download/elastix-5.0.0-manual.pdf
6.  ElastixFromMatlab (a MATLAB® wrapper around elastix)
  by: CNRS,France and Riverside Research, USA https://sourcesup.renater.fr/www/elxfrommatlab/
  - It is included in the Dependency folder, so no need to download it.
  - In case of you had to redownload the elastix toolbox as mentioned above, make sure to change the paths in "elxTestDefaultConfiguration.m" script.
7.  loadImage3
  by: Dr. Mattijs Elschot from the MR center at the Norwegian University of Science and Technology (NTNU), Trondheim, Norway.
  Dr. Elschot allowed the function useage and upload. 
8. SegmentationQualityControl method
  by MR Cancer Group at NTNU Trondheim, Norway
  - It is included in the code, no need to install it, make sure it works and that you have all of its dependencies as descriped in:
  https://github.com/ntnu-mr-cancer/SegmentationQualityControl 

# Contact us
Feel free to contact us:
mohammed.sunoqrot@ntnu.no

