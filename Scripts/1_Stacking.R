#We first load the required R packages:
library(data.table)
library(celestial)
library(Rfits)
library(ProFound)
library(doParallel)
library(foreach)
library(imager)
library(ProPane)

registerDoParallel(cores=8) #parallel execution of R code on machines with multiple cores (8, in my case)

#For the purpose of this work, we only use the overlapped data between HSC-PDR3 and GAMA. We first reconstruct 
#a PSF per GAMA region and band:

region = 'G02' #then: G09, G12, G15
band = 'g' #then: r,i,Z,Y

#We select the stars depending of its magnitude for each of the 4 PSF regions. GAIA_region.csv is the GAIA
#star catalogue from each GAMA region, located at GAIA_catalogues/
GAIA = fread(paste0('/Path2GAIACatalogue/GAIA_',region,'.csv'))
GAIA_bright = GAIA[phot_g_mean_mag < 8,] 
GAIA_mod = GAIA[phot_g_mean_mag > 11 & phot_g_mean_mag < 11.5,]
GAIA_faint = GAIA[phot_g_mean_mag > 14 & phot_g_mean_mag < 14.1,]
GAIA_faintest = GAIA[phot_g_mean_mag > 18 & phot_g_mean_mag < 18.02,]

#We read HSC data previoulsy downloaded (select filter) and create a .csv file w/info of each HSC fits file in
#the overlapped region between GAMA and HSC:
suppressMessages({
temp_scan = Rfits_key_scan(dirlist = '/Path2HSCData/hsc-release.mtk.nao.ac.jp/archive/filetree/pdr3_wide/deepCoadd/HSC-G', cores=8, get_all=TRUE, ext=2)
})
write.csv(temp_scan, paste0('/Path2HSCFileInfo/file_info_',band,'.csv')) 

temp_scan = read.csv(paste0('/Path2HSCFileInfo/file_info_',band,'.csv')) 

#We make the match between the GAIA catalogue with info of the stars and the HSC fits files, where the matching radius is 1 deg.
#refID gives the row position from GAIA catalogue and compareID gives the corresponding best matching row position in HSC:
good_match_bright = coordmatch(GAIA_bright[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')
good_match_mod = coordmatch(GAIA_mod[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')
good_match_faint = coordmatch(GAIA_faint[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')
good_match_faintest = coordmatch(GAIA_faintest[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


                                       ###########OUTER STACK###########
dummy = foreach(ID = good_match_bright$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_bright$ID[ID,1],"full"], ext=2) #We read the corresponding HSC image in each pf the matches
  good_star_cut = good_star[GAIA_bright[ID,ra], GAIA_bright[ID,dec], box = 4001, type='coord'] #We center the HSC image in the GAIA star
  png(paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  #We create the fits files with the unmasked background sources that we'll stack later
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits')) 
  
  #We mask the background sources only for the normalisation step. I work with a smaller cut that contains the normalization annulus (400-410 pix)     
  good_star_cut = good_star[GAIA_bright[ID,ra], GAIA_bright[ID,dec], box = 1001, type='coord']
  
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
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=7) > 0.65
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=3) 
  }
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0 | is.na(good_star_cut))) / prod(dim(good_star_cut_mask)) #Creat keyvalue MASK with proportion of masked pixels
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  #Save masked fits file
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  Rfits_write_image(good_star_cut_mask, paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'), create_file=FALSE, create_ext=TRUE)
  
  good_star_cut$imDat[good_star_cut_mask > 0] = NA #Puts NA where we have the mask
  
  png(paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}

#We read the info from the fits files we just created
bright_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

#We expand the grid that has only x,y positions from each pixel of the image. We add a third column with the radius measured from
#the centre of the image to each pixel:
bright_rad = (bright_scan[1,"dim_1"][[1]] -1L) / 2
sel_grid = expand.grid(-bright_rad:bright_rad, -bright_rad:bright_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 400 & sel_grid[,3] < 410) #this is the normalization annulus! We normalize the star images before stacking

#We make a list with all the bright stars that we are going to stack:
bright_list = foreach(ID = good_match_bright$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  
  temp_im = Rfits_read(paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Bright_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.8){ #we discard the ones that have more than 80% of the pixels masked
    temp_mask[[1]]$imDat[which(temp_mask[[2]]$imDat > 0)] = NA
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    temp_im[[1]]$imDat = temp_im[[1]]$imDat/scale #normalisation step
    temp_im[[1]]$imDat[is.nan(temp_im[[1]]$imDat)] = NA
    return(temp_im[[1]]$imDat)
  }else{
    return(stop('Skipping', ID))
  }
}

#Median stack with the unmasked images:
stack_bright = propaneStackFlatFunc(image_list = bright_list, imager_func = parmed)

stack_bright_sd = propaneStackFlatFunc(image_list = bright_list, imager_func = parsd)

rm(bright_list)

stack_bright$stdev = stack_bright_sd$image

Rfits_write(stack_bright, filename = paste0('/Path2StarFits/Bright_Star_Stack_',region,'_',band,'.fits'))


                                       ###########MIDDLE STACK###########

dummy = foreach(ID = good_match_mod$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_mod$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_mod[ID,ra], GAIA_mod[ID,dec], box = 1001, type='coord']
  
  png(paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
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
  
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0 | is.na(good_star_cut))) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  Rfits_write_image(good_star_cut_mask, paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'), create_file=FALSE, create_ext=TRUE)
  
  good_star_cut$imDat[good_star_cut_mask > 0] = NA
  
  png(paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}

mod_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

mod_rad = (mod_scan[1,"dim_1"][[1]] -1L) / 2
sel_grid = expand.grid(-mod_rad:mod_rad, -mod_rad:mod_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 200 & sel_grid[,3] < 220)

mod_list =  foreach(ID = good_match_mod$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_im = Rfits_read(paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Mod_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.5){
    temp_mask[[1]]$imDat[which(temp_mask[[2]]$imDat > 0)] = NA
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    temp_im[[1]]$imDat = temp_im[[1]]$imDat/scale
    temp_im[[1]]$imDat[is.nan(temp_im[[1]]$imDat)] = NA
    return(temp_im[[1]]$imDat)
  }else{
    return(stop('Skipping', ID))
  }
}

stack_mod = propaneStackFlatFunc(image_list = mod_list, imager_func = parmed)

stack_mod_sd = propaneStackFlatFunc(image_list =  mod_list, imager_func = parsd)

rm(mod_list)
stack_mod$stdev = stack_mod_sd$image

Rfits_write(stack_mod, filename = paste0('/Path2StarFits/Mod_Star_Stack_',region,'_',band,'.fits'))


                                       ###########INNER STACK###########

dummy = foreach(ID = good_match_faint$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_faint$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_faint[ID,ra], GAIA_faint[ID,dec], box = 501, type='coord']
  
  png(paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
  if(band == 'g' | band == 'r'| band == 'i'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=1) > 0.2
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=3)
  }
  
  if(band == 'Y'| band == 'Z'){
    good_star_cut_diff = profoundImDiff(good_star_cut$imDat, sigma=2) > 0.9
    good_star_cut_mask = profoundDilate(good_star_cut_diff, size=2)
  }
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0 | is.na(good_star_cut))) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  Rfits_write_image(good_star_cut_mask, paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'), create_file=FALSE, create_ext=TRUE)
  
  good_star_cut$imDat[good_star_cut_mask > 0] = NA
  
  png(paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  
  return(NULL)
}


faint_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

faint_rad = (faint_scan[1,"dim_1"][[1]] -1L) / 2 
sel_grid = expand.grid(-faint_rad:faint_rad, -faint_rad:faint_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 100 & sel_grid[,3] < 140)

faint_list = foreach(ID = good_match_faint$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_im = Rfits_read(paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Faint_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.02){
    temp_mask[[1]]$imDat[which(temp_mask[[2]]$imDat > 0)] = NA
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    temp_im[[1]]$imDat = temp_im[[1]]$imDat/scale
    temp_im[[1]]$imDat[is.nan(temp_im[[1]]$imDat)] = NA
    return(temp_im[[1]]$imDat)
  }else{
    return(stop('Skipping', ID))
  }
}

stack_faint = propaneStackFlatFunc(image_list = faint_list, imager_func = parmed)

stack_faint_sd = propaneStackFlatFunc(image_list =  faint_list, imager_func = parsd)

rm(faint_list)

stack_faint$stdev = stack_faint_sd$image

Rfits_write(stack_faint, filename = paste0('/Path2StarFits/Faint_Star_Stack_',region,'_',band,'.fits'))


                                       ###########COREE STACK###########

dummy = foreach(ID = good_match_faintest$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_faintest$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_faintest[ID,ra], GAIA_faintest[ID,dec], box = 201, type='coord']
  
  png(paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
  segim = profoundProFound(good_star_cut, box=50)$segim #Run profound and select segmentation maps
 
  good_star_cut_mask = segim > 0 | is.na(good_star_cut) #Mask where we have positive segim or NA
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0)) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  Rfits_write_image(good_star_cut_mask, paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'), create_file=FALSE, create_ext=TRUE)
  
  good_star_cut$imDat[good_star_cut_mask > 0] = NA
  
  png(paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}


faintest_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

faintest_rad = (faintest_scan[1,"dim_1"][[1]] -1L) / 2 
sel_grid = expand.grid(-faintest_rad:faintest_rad, -faintest_rad:faintest_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 50 & sel_grid[,3] < 100)


faintest_list = foreach(ID = good_match_faintest$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_im = Rfits_read(paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Faintest_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.2){
    temp_mask[[1]]$imDat[which(temp_mask[[2]]$imDat > 0)] = NA
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    if(which.max(temp_im[[1]]$imDat) != (faintest_rad*(faintest_rad*2 + 1) + faintest_rad + 1)){
      #We only select the faintest stars that are exactly centred! This is, the stars that have the brightest pixel in the center
      return(stop('Skipping', ID))
    }
    temp_im[[1]]$imDat = temp_im[[1]]$imDat/scale
    temp_im[[1]]$imDat[is.nan(temp_im[[1]]$imDat)] = NA
    return(temp_im[[1]]$imDat)
  }else{
    return(stop('Skipping', ID))
  }
}

stack_faintest = propaneStackFlatFunc(image_list = faintest_list, imager_func = parmed)

stack_faintest_sd = propaneStackFlatFunc(image_list =  faintest_list, imager_func = parsd)

rm(faintest_list)

stack_faintest$stdev = stack_faintest_sd$image

Rfits_write(stack_faintest, filename = paste0('/Path2StarFits/Faintest_Star_Stack_',region,'_',band,'.fits'))
