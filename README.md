# HSC-SSP PDR3 PSFs

## PSF fits files

The HSC-SSP PDR3 PSF 2-D models are available as FITS files at PSF_fits/.
To download a single file, start by clicking the **Go to file** button at the top of the repository contents. 
This will pull up a page that lists all of the files in the HSC_PSFs repository. 
Select the folder PSF_fits and click on the file you wish to download to open the individual file. From here, click the **Raw** button at the top of the file, this will automatically save the file.

To download the complete repository, click the green **Code** button at the top right of the repo contents and select Download Zip. This will download the entire GitHub repository as a compressed zipped folder. 

##  Scripts for PSF reconstruction
The reconstruction of the PSF is divided into 3 scripts available at Scrptis/: 1_Stacking.R, 2_Combining.R, and 3_FinalPSF.R. 1_Stacking.R selects the coordinates of stars from the GAIA catalogue (available as CSV files at GAIA_catalogues/), finds the HSC image of each star, and performs the stack of the images. This is done for each region of the PSF (outer, middle, inner, and core), per GAMA region (G02, G09, G12, and G15), and per HSC band (g, r, i, Z, Y). The resulting stack is saved as a FITS file. The next step is 2_Combining.R, which combines the outer, middle, inner, and core parts to reconstruct a PSF pero GAMA region and HSC band. the resulting PSF is also saved as a FITS file. Finally, 3_FinalPSF.R stacks the 4 HSC PSFs per band to create the final HSC-SSP PDR3 PSFs.


##  Install R 
All scripts are written in the open-source **R** language. To download a recent version of **R**, go to the **R** project page <https://cloud.r-project.org/> and select the binary distribution depending on your operating system. Then:

#### Mac
Click on *R-4.3.1-arm64.pkg* (or newest version) to download R. Once the file download is complete, click to open the installer. Click Continue and proceed through the installer and follow the instructions.

#### Windows
Select the subdirectory **base** if you are installing **R** for the first time, then **Download R-4.3.1 for Windows**. The distribution is distributed as an installer *R-4.3.1-win.exe*, run this file for a Windows-style installer.

#### Linux
Select the Linux distribution. The exact installation procedure will vary depending on the Linux system you use. CRAN guides the process by grouping each set of source files with documentation or README files that explain how to install on your system.

##  Install R-studio
Once the  **R** installer has finished, the recommendation is to use the **R** integrated development environment (IDE) **R-studio**. Its interface is well organized so that the user can clearly view graphs, data tables, R code, and output all at the same time. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/>.


##  R packages:
The two required packages that are not available to download from CRAN are **ProFound** and **ProPane**:
#### ProFound
#### ProPane








