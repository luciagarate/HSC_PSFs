# HSC-SSP PDR3 PSFs

## PSF fits files

The HSC-SSP PDR3 PSF 2-D models are available as FITS files at PSF_fits/.
To download a single file, start by clicking the **Go to file** button at the top of the repository contents. 
This will pull up a page that lists all of the files in the HSC_PSFs repository. 
Select the folder PSF_fits and click on the file you wish to download to open the individual file. From here, click the **Raw** button at the top of the file, this will automatically save the file.

To download the complete repository, use the command line *git-lfs clone <repository.>* Git LFS is a command line extension and specification for managing large files with Git, visit <https://github.com/git-lfs/git-lfs> to install it. This will download the entire GitHub repository as a compressed zipped folder. 

The PSF FITS files are compressed FITS files of 4001x4001 pix<sup>2</sup>, where the second extension contains the image data.

##  Scripts for PSF reconstruction
The reconstruction of the PSF is divided into 3 scripts available at Scripts/: "1_StarSelection.R", "2_Stacking.R", and "3_Combining.R". "1_StarSelection.R" selects the coordinates of stars from the GAIA catalogue (available as CSV files at GAIA_catalogues/), finds the HSC image of each star, and creates cutouts. This is done for each region of the PSF (outer, middle, inner, and core) and per HSC band (g, r, i, Z, Y). The next step is "2_Stacking.R", which performs an stack of all the selected stars per PSF region and band.  In the last script, "3_Combining.R", we construct an HSC-SSP PDR3PSF per HSC band by combining the outer + middle + inner + core regions. The resulting stack is saved as a FITS file.


##  Install R 
All scripts are written in the open-source **R** language. To download a recent version of **R**, go to the **R** project page <https://cloud.r-project.org/> and select the binary distribution depending on your operating system. Then:

#### Mac
Click on *R-4.3.1-arm64.pkg* (or newest version) to download R. Once the file download is complete, click to open the installer. Click **Continue** and proceed through the installer.

#### Windows
Select the subdirectory **base** if you are installing **R** for the first time, then **Download R-4.3.1 for Windows**. The distribution is distributed as an installer *R-4.3.1-win.exe*, run this file for a Windows-style installer.

#### Linux
Select the **Linux distribution**. The exact installation procedure will vary depending on the Linux system you use. CRAN guides the process by grouping each set of source files with documentation or README files that explain how to install on your system.

##  Install R-studio
Once the  **R** installer has finished, the recommendation is to use the **R** integrated development environment (IDE) **R-studio**. Its interface is well organized so that the user can clearly view graphs, data tables, R code, and output all at the same time. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/>.


##  R packages:
The two required packages that are not available to download from CRAN are **ProFound** and **ProPane**:
#### ProFound

Source installation from GitHub:

```R
install.packages('devtools')
install_github("asgr/ProFound")
library(ProFound)
```
Or visit <https://github.com/asgr/ProFound/>

#### ProPane

Source installation from GitHub:

```R
install.packages('devtools')
install_github("asgr/ProPane")
library(ProFound)
```
Or visit <https://github.com/asgr/ProPane/>








