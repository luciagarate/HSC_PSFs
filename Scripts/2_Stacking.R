library(data.table)
library(celestial)
library(Rfits)
library(ProFound)
library(doParallel)
library(foreach)


band = 'g'

#Same as in "1_StarSelection.R":
temp_scan = read.csv(paste0('/Path2HSCFileInfo/file_info_',band,'.csv')) #tris was created in "1_StarSelection.R"
GAIA_outer = fread(paste0('/Path2GAIACatalogue/GAIA_OuterStars_HSCWIDE.csv'))
good_match_outer = coordmatch(GAIA_outer[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


#We expand the grid that has only x,y positions from each pixel of the image. We add a third column with the radius measured from
#the centre of the image to each pixel:
outer_rad = 2000 #(Size of one outer cutout - 1)/2
sel_grid = expand.grid(-outer_rad:outer_rad, -outer_rad:outer_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 400 & sel_grid[,3] < 410) #this is the normalization annulus! We normalize the star images before stacking

#We save in a folder all the masked normalised bright stars that we are going to stack, so in this loop we first discard the 
#images that have more than 80% of the pixels masked and we then normalise them:
foreach(ID = good_match_outer$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_mask = Rfits_read(paste0('/Path2StarFits/Outer_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits')) #Read the masked fits file
  if(temp_mask[[1]]$keyvalues$MASK < 0.8){ #we discard the ones that have more than 80% of the pixels masked
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    temp_mask[[1]]$imDat = temp_mask[[1]]$imDat/scale #normalisation step
    temp_mask[[1]]$imDat[is.nan(temp_mask[[1]]$imDat)] = NA
    Rfits_write_image(temp_mask[,]$imDat,paste0('/Path2StarFits/Outer_Stars_',band,'/2Stack_masked/GAIA_',formatC(ID,width=3,flag=0),'_norm_mask.fits'))
  }else{
    return(stop('Skipping', ID))
  }
}

                                

#Median stack with the masked images, in this case we select cores=4 and chunk=501, so it divides the 4001x4001 stack in 64 regions
stack_outer = propaneStackWarpMed(dirlist='/Path2StarFits/Outer_Stars_',band,'/2Stack_masked/', useCUTLO = FALSE, keyvalues_out = NULL, chunk=501, cores = 4)
Rfits_write(stack_outer, filename = paste0('/Path2Stacks/Outer_Star_Stack_',band,'.fits'))

                                            ###########MIDDLE STACK###########


GAIA = fread(paste0('/Path2GAIACatalogue/GAIA_all_stars_HSC_GAMA.csv'))
GAIA_mid = GAIA[phot_g_mean_mag > 11 & phot_g_mean_mag < 11.5,]
good_match_mid = coordmatch(GAIA_mid[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


mid_rad = (1001-1)/2 #(Size of one middle cutout - 1)/2
sel_grid = expand.grid(-mid_rad:mid_rad, -mid_rad:mid_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 200 & sel_grid[,3] < 220)

foreach(ID = good_match_mid$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_mask = Rfits_read(paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.5){
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    temp_mask[[1]]$imDat = temp_mask[[1]]$imDat/scale
    temp_mask[[1]]$imDat[is.nan(temp_mask[[1]]$imDat)] = NA
    Rfits_write_image(temp_mask[,]$imDat,paste0('/Path2StarFits/Mid_Stars_',band,'/2Stack_masked/GAIA_',formatC(ID,width=3,flag=0),'_norm_mask.fits'))
  }else{
    return(stop('Skipping', ID))
  }
}

stack_mid = propaneStackWarpMed(dirlist='/Path2StarFits/Mid_Stars_',band,'/2Stack_masked/', useCUTLO = FALSE, keyvalues_out = NULL, cores = 4)
Rfits_write(stack_outer, filename = paste0('/Path2Stacks/Mid_Star_Stack_',band,'.fits'))



                                      ###########INNER STACK###########
GAIA_inner = GAIA[phot_g_mean_mag > 14 & phot_g_mean_mag < 14.1,]
good_match_inner = coordmatch(GAIA_inner[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')

inner_rad = (501-1)/2 #(Size of one inner cutout - 1)/2 
sel_grid = expand.grid(-inner_rad:inner_rad, -inner_rad:inner_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 100 & sel_grid[,3] < 140)

inner_list = foreach(ID = good_match_inner$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_mask = Rfits_read(paste0('/Path2StarFits/Inner_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.02){
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    temp_mask[[1]]$imDat = temp_mask[[1]]$imDat/scale #normalisation
    temp_mask[[1]]$imDat[is.nan(temp_mask[[1]]$imDat)] = NA
    Rfits_write_image(temp_mask[,]$imDat,paste0('/Path2StarFits/Inner_Stars_',band,'/2Stack_masked/GAIA_',formatC(ID,width=3,flag=0),'_norm_mask.fits'))
  }else{
    return(stop('Skipping', ID))
  }
}

stack_inner = propaneStackWarpMed(dirlist='/Path2StarFits/Inner_Stars_',band,'/2Stack_masked/', useCUTLO = FALSE, keyvalues_out = NULL, cores=4)
Rfits_write(stack_inner, filename = paste0('/Path2Stacks/Inner_Star_Stack_',band,'.fits'))

                                      ###########CORE STACK###########


GAIA_core = GAIA[phot_g_mean_mag > 18 & phot_g_mean_mag < 18.02,]
good_match_core = coordmatch(GAIA_core[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')



core_rad = (200-1)/2 #(Size of one core cutout - 1)/2
sel_grid = expand.grid(-core_rad:core_rad, -core_rad:core_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 0 & sel_grid[,3] < 100)


foreach(ID = good_match_core$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_im = Rfits_read(paste0('/Path2StarFits/Core_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Core_Stars_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.2){
    scale = median(temp_mask[[1]]$imDat[sel_pix], na.rm = TRUE)
    if(scale < 0){
      return(stop('Skipping', ID))
    }
    if(which.max(temp_im[[1]]$imDat) != (core_rad*(core_rad*2 + 1) + core_rad + 1)){
      #We only select the faintest stars that are exactly centred! This is, the stars that have the brightest pixel in the center
      return(stop('Skipping', ID))
    }
    temp_im[[1]]$imDat = temp_im[[1]]$imDat/scale
    temp_im[[1]]$imDat[is.nan(temp_im[[1]]$imDat)] = NA
    Rfits_write_image(temp_mask[,]$imDat,paste0('/Path2StarFits/Core_Stars_',band,'/2Stack_masked/GAIA_',formatC(ID,width=3,flag=0),'_norm_mask.fits'))
  }else{
    return(stop('Skipping', ID))
  }
}

stack_core = propaneStackWarpMed(dirlist='/Path2StarFits/Core_Stars_',band,'/2Stack_masked/', useCUTLO = FALSE, keyvalues_out = NULL, cores=4)
Rfits_write(stack_core, filename = paste0('/Path2Stacks/Core_Star_Stack_',band,'.fits'))