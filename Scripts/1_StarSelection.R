#We first load the required R packages:
library(data.table)
library(celestial)
library(Rfits)
library(ProFound)
library(doParallel)
library(foreach)
library(imager) #parmed function in propaneStackFlatFunc
library(ProPane)

registerDoParallel(cores=8) #parallel execution of R code on machines with multiple cores (8, in my case)

band = 'g' #then: r,i,Z,Y

#For the purpose of this work, we only use all HSC Wide data for the outer PSF region and the overlapped data between HSC-PDR3 and GAMA for
#the inner, middle and core regions.

                                              ###########OUTER STACK###########
GAIA_outer = fread(paste0('/Path2GAIACatalogue/GAIA_OuterStars_HSCWIDE.csv')) #This file contains information about all the stars with mag<8
#in all HSC Wide, which are 1739 stars


#We read HSC data previoulsy downloaded (select filter) and create a .csv file w/info of each HSC fits file in
#the overlapped region between GAMA and HSC:
suppressMessages({
  temp_scan = Rfits_key_scan(dirlist = '/Path2HSCData/hsc-release.mtk.nao.ac.jp/archive/filetree/pdr3_wide/deepCoadd/HSC-G', cores=8, get_all=TRUE, ext=2)
})
write.csv(temp_scan, paste0('/Path2HSCFileInfo/file_info_',band,'.csv')) 

temp_scan = read.csv(paste0('/Path2HSCFileInfo/file_info_',band,'.csv')) 


#We make the match between the GAIA catalogue with info of the stars and the HSC fits files, where the matching radius is 1 deg.
#refID gives the row position from GAIA catalogue and compareID gives the corresponding best matching row position in HSC:
good_match_outer = coordmatch(GAIA_outer[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


dummy = foreach(ID = good_match_outer$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_outer$ID[ID,1],"full"], ext=2) #We read the corresponding HSC image in each of the matches
  good_star_cut = good_star[GAIA_outer[ID,ra], GAIA_outer[ID,dec], box = 4001, type='coord'] #We center the HSC image in the GAIA star and make cutout
  
  png(paste0('/Path2StarFits/Outer_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Outer_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits')) 
  
  #We mask the background sources:
  #profoundImDiffunction creates an image which is the original minus a smoothed version.
  #Sigma establishes the standard deviation of the blur. profoundDilate dilates the seg maps with the parameter size. By only selecting the 
  #pixels with profoundImDiff > value, we decide how aggressive we want to be with the masks.
  #g,r,i bands:
  if(band == 'g' | band == 'r'| band == 'i'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=3) > 0.2
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=21) 
  }
  
  #Z,Y bands:
  if(band == 'Z' | band == 'Y'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=25) > 1
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=13)
  }
  
  outer_rad = (4001 -1L) / 2 
  sel_grid = expand.grid(-outer_rad:outer_rad, -outer_rad:outer_rad)
  sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
  sel_pix = which(sel_grid[,3] <  200)
  
  if(band == 'r'){
    bright_rad = (4001 -1L) / 2 
    sel_grid = expand.grid(-bright_rad:bright_rad, -bright_rad:bright_rad)
    sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
    sel_pix = which(sel_grid[,3] <  340)
  }
  
  
  original = good_star_cut$imDat
  good_star_cut$imDat[good_star_cut_mask > 0] = NA #Puts NA where we have the mask
  good_star_cut$imDat[sel_pix] = original[sel_pix] #Unmask pixels in a R<200 pix
  
  
  sel_pix_horiz = which(sel_grid[,2] <  40 & sel_grid[,2] > -40) #Unmask the horizontal spikes in all the bands
  sel_pix_vert = which(sel_grid[,1] <  40 & sel_grid[,1] > -40) #Unmask the vertical spikes in the Y-band
  
  good_star_cut$imDat[sel_pix_horiz] = original[sel_pix_horiz]
  
  if(band == 'Y'){
    good_star_cut$imDat[sel_pix_vert] = original[sel_pix_vert]
  }
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0 | is.na(good_star_cut_mask))) / prod(dim(good_star_cut_mask)) #Creat keyvalue MASK with proportion of masked pixels
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  #Save masked fits file
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Outer_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))

  png(paste0('/Path2StarFits/Outer_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}



                                          ###########MIDDLE STACK###########

#As said before, for the middle, inner and core stacks the stars that are in the overlapped region between HSC and GAMA are enough,
#so GAIA_all_stars_HSC_GAMA.csv has the information of all the GAIA stars that are in these regions - we have to select the ones that we want depending on the magnitude range:
#GAIA_all_stars_HSC_GAMA.csv is located in GAIA_catalogues/
GAIA = fread(paste0('/Path2GAIACatalogue/GAIA_all_stars_HSC_GAMA.csv'))
GAIA_mid = GAIA[phot_g_mean_mag > 11 & phot_g_mean_mag < 11.5,]
good_match_mid = coordmatch(GAIA_mid[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


dummy = foreach(ID = good_match_mid$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_mid$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_mid[ID,ra], GAIA_mid[ID,dec], box = 1001, type='coord']
  
  png(paste0('/Path2StarFits/Mid_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Mid_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
  if(band == 'g' | band == 'r'| band == 'i'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=2) > 0.2
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=5)
  }
  if(band == 'Z'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=7) > 0.65
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=3) 
  }
  if(band == 'Y'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=13) > 0.9
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=9)
  }
  
  
  mid_rad = (1001 -1L) / 2 
  sel_grid = expand.grid(-mid_rad:mid_rad, -mid_rad:mid_rad)
  sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
  sel_pix = which(sel_grid[,3] <  150) #Select pixels inside a R = 150 pix
  
  original = good_star_cut$imDat
  good_star_cut$imDat[good_star_cut_mask > 0] = NA #Puts NA where we have the mask
  good_star_cut$imDat[sel_pix] = original[sel_pix] #Unmask pixels in a R<150 pix
  
  
  sel_pix_horiz = which(sel_grid[,2] <  40 & sel_grid[,2] > -40) #Unmask the horizontal spikes in all the bands
  sel_pix_vert = which(sel_grid[,1] <  40 & sel_grid[,1] > -40) #Unmask the vertical spikes in the Y-band
  
  good_star_cut$imDat[sel_pix_horiz] = original[sel_pix_horiz]
  
  if(band == 'Y'){
    good_star_cut$imDat[sel_pix_vert] = original[sel_pix_vert]
  }
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0 | is.na(good_star_cut_mask))) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Mid_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  
  
  png(paste0('/Path2StarFits/Mid_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}

                                       ###########INNER STACK###########

GAIA_inner = GAIA[phot_g_mean_mag > 14 & phot_g_mean_mag < 14.1,]
good_match_inner = coordmatch(GAIA_inner[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


dummy = foreach(ID = good_match_inner$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_inner$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_inner[ID,ra], GAIA_inner[ID,dec], box = 501, type='coord']
  
  png(paste0('/Path2StarFits/Inner_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Inner_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
  if(band == 'g' | band == 'r'| band == 'i'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=1) > 0.2
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=3)
  }
  
  if(band == 'Y'| band == 'Z'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=2) > 0.9
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=2)
  }
  
  
   original = good_star_cut$imDat
   good_star_cut$imDat[good_star_cut_mask > 0] = NA
   
   #Unmask vertical spike in Y-band
   if(band == 'Y'){
       inner_rad = (501 -1L) / 2
       sel_grid = expand.grid(-inner_rad:inner_rad, -inner_rad:inner_rad)
       sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)

       sel_pix_horiz = which(sel_grid[,2] <  12 & sel_grid[,2] > -12)
       sel_pix_vert = which(sel_grid[,1] <  12 & sel_grid[,1] > -12)

       good_star_cut$imDat[sel_pix_vert] = original[sel_pix_vert]
       good_star_cut$imDat[sel_pix_horiz] = original[sel_pix_horiz]

   }
   
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0 | is.na(good_star_cut_mask))) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Inner_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  
  
  png(paste0('/Path2StarFits/Inner_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  
  return(NULL)
}


                                       ###########CORE STACK###########

GAIA_core = GAIA[phot_g_mean_mag > 18 & phot_g_mean_mag < 18.02,]
good_match_core = coordmatch(GAIA_core[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


dummy = foreach(ID = good_match_core$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_core$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_core[ID,ra], GAIA_core[ID,dec], box = 201, type='coord']
  
  png(paste0('/Path2StarFits/Core_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Core_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
  segim = profoundProFound(good_star_cut, box=50)$segim #Run profound and select segmentation maps
  
  seg_cen = segim[101,101] #Select the central segment
  segim[segim == seg_cen] = 0L 
  
  good_star_cut_mask = segim > 0 | is.na(good_star_cut) 
  
  good_star_cut$imDat[good_star_cut_mask > 0] = NA
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0)) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  

  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Core_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  
  png(paste0('/Path2StarFits/Core_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}

