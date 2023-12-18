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

#For the purpose of this work, we only use the overlapped data between HSC-PDR3 and GAMA. We first reconstruct 
#a PSF per GAMA region and band:

region = 'G02' #then: G09, G12, G15
band = 'g' #then: r,i,Z,Y

#We select the stars depending of its magnitude for each of the 4 PSF regions. GAIA_region.csv is the GAIA
#star catalogue from each GAMA region, located at GAIA_catalogues/
GAIA = fread(paste0('/Path2GAIACatalogue/GAIA_',region,'.csv'))
GAIA_outer = GAIA[phot_g_mean_mag < 8,] 
GAIA_mid = GAIA[phot_g_mean_mag > 11 & phot_g_mean_mag < 11.5,]
GAIA_inner = GAIA[phot_g_mean_mag > 14 & phot_g_mean_mag < 14.1,]
GAIA_core = GAIA[phot_g_mean_mag > 18 & phot_g_mean_mag < 18.02,]

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
good_match_mid = coordmatch(GAIA_mid[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')
good_match_inner = coordmatch(GAIA_inner[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')
good_match_core = coordmatch(GAIA_core[,list(ra,dec)], temp_scan[,c("centre_RA", "centre_Dec")], rad=1, radunit = 'deg')


                                       ###########OUTER STACK###########
dummy = foreach(ID = good_match_outer$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_outer$ID[ID,1],"full"], ext=2) #We read the corresponding HSC image in each pf the matches
  good_star_cut = good_star[GAIA_outer[ID,ra], GAIA_outer[ID,dec], box = 4001, type='coord'] #We center the HSC image in the GAIA star and make cutout
  png(paste0('/Path2StarFits/Outer_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  #We create the fits files with the unmasked background sources that we'll stack later
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Outer_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits')) 
  
  #We mask the background sources only for the normalisation step:
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
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Outer_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  
  
  png(paste0('/Path2StarFits/Outer_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}


list_IDs = good_match_bright$bestmatch[,1]
#Manually iscar bad star images of the outer stack that could introduce error to the outer stack:

if(region == 'G12'){
  index_delete = which(list_IDs == 1 | list_IDs == 2 | list_IDs == 3| list_IDs == 10 | list_IDs == 23 | list_IDs == 24 | list_IDs == 44 | list_IDs == 49 | list_IDs == 51 | list_IDs == 52 | list_IDs == 57 | list_IDs == 61)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G12' & band == 'r'){
  index_delete = which(list_IDs == 5 | list_IDs == 9 | list_IDs == 18 | list_IDs == 19 | list_IDs == 21 | list_IDs == 32 | list_IDs == 40 | list_IDs == 56 | list_IDs == 64 | list_IDs == 65)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G12' & band == 'Y'){
  index_delete = which(list_IDs == 5 | list_IDs == 9 | list_IDs == 18 | list_IDs == 19 | list_IDs == 32 | list_IDs == 34 | list_IDs == 40 | list_IDs == 56 | list_IDs == 68)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G12' & band == 'Z'){
  index_delete = which(list_IDs == 5 | list_IDs == 9 | list_IDs == 18 | list_IDs == 19 | list_IDs == 21 | list_IDs == 32 | list_IDs == 40 | list_IDs == 48 | list_IDs == 56 | list_IDs == 58 | list_IDs == 65 | list_IDs == 68)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G02'){
  index_delete = which(list_IDs == 1 | list_IDs == 12 | list_IDs == 25)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G02'  & band == 'Y'){
  index_delete = which(list_IDs == 5)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G02'  & band == 'Z'){
  index_delete = which(list_IDs == 5 | list_IDs == 12 | list_IDs == 25)
  list_IDs = list_IDs[-rbind(index_delete)]
}


if(region == 'G09'){
  index_delete = which(list_IDs == 1 | list_IDs == 4 | list_IDs == 5 | list_IDs == 21 | list_IDs == 25 | list_IDs == 60 | list_IDs == 68 | list_IDs == 79 | list_IDs == 90)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G09' & band == 'Y'){
  index_delete = which(list_IDs == 11 | list_IDs == 23 | list_IDs == 93)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G09' & band == 'r'){
  index_delete = which(list_IDs == 11 | list_IDs == 54 | list_IDs == 70 | list_IDs == 80 )
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G09' & band == 'Z'){
  index_delete = which(list_IDs == 6 | list_IDs == 11 | list_IDs == 17 | list_IDs == 22 | list_IDs == 23 | list_IDs == 29| list_IDs == 55 | list_IDs == 57 | list_IDs == 65 | list_IDs == 93)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G15'){
  index_delete = which(list_IDs == 21 | list_IDs == 23  | list_IDs == 26 | list_IDs == 32 | list_IDs == 38 | list_IDs == 47 | list_IDs == 54 | list_IDs == 67)
  list_IDs = list_IDs[-rbind(index_delete)]
}


if(region == 'G15' & band == 'r'){
  index_delete = which(list_IDs == 2 | list_IDs == 9 | list_IDs == 14 |list_IDs == 19 | list_IDs == 21 | list_IDs == 23  | list_IDs == 26 |list_IDs == 27 |list_IDs == 28 | list_IDs == 31 |list_IDs == 32 | list_IDs == 33 |list_IDs == 36 |list_IDs == 38 | list_IDs == 39 |list_IDs == 42  | list_IDs == 45  |list_IDs == 47 | list_IDs == 54 | list_IDs == 57  |list_IDs == 67)
  list_IDs = list_IDs[-rbind(index_delete)]
}


if(region == 'G15' & band == 'Y'){
  index_delete = which(list_IDs == 2 | list_IDs == 14 |list_IDs == 17 | list_IDs == 19 |list_IDs == 27 | list_IDs == 28  | list_IDs == 30 |list_IDs == 33 |list_IDs == 34 | list_IDs == 39 |list_IDs == 45 | list_IDs == 41)
  list_IDs = list_IDs[-rbind(index_delete)]
}

if(region == 'G15' & band == 'Z'){
  index_delete = which(list_IDs == 1 | list_IDs == 2 |list_IDs == 9 | list_IDs == 14 |list_IDs == 17 | list_IDs == 57)
  list_IDs = list_IDs[-rbind(index_delete)]
}



#We read the info from the fits files we just created
outer_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Outer_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

#We expand the grid that has only x,y positions from each pixel of the image. We add a third column with the radius measured from
#the centre of the image to each pixel:
outer_rad = (outer_scan[1,"dim_1"][[1]] -1L) / 2
sel_grid = expand.grid(-outer_rad:outer_rad, -outer_rad:outer_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 400 & sel_grid[,3] < 410) #this is the normalization annulus! We normalize the star images before stacking

#We make a list with all the bright stars that we are going to stack:
outer_list = foreach(ID = list_IDs, .errorhandling = 'remove')%dopar%{
  
  temp_im = Rfits_read(paste0('/Path2StarFits/Outer_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Outer_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.8){ #we discard the ones that have more than 80% of the pixels masked
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
stack_outer = propaneStackFlatFunc(image_list = outer_list, imager_func = parmed)

stack_outer_sd = propaneStackFlatFunc(image_list = outer_list, imager_func = parsd)

rm(outer_list)

stack_outer$stdev = stack_outer_sd$image

Rfits_write(stack_outer, filename = paste0('/Path2StarFits/Outer_Star_Stack_',region,'_',band,'.fits'))


                                       ###########MIDDLE STACK###########

dummy = foreach(ID = good_match_mid$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_mid$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_mid[ID,ra], GAIA_mid[ID,dec], box = 1001, type='coord']
  
  png(paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
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
  
  
 original = good_star_cut$imDat
 good_star_cut$imDat[good_star_cut_mask > 0] = NA

 #Unmask vertical spike in Y-band
 if(band == 'Y'){
 mid_rad = (1001 -1L) / 2
 sel_grid = expand.grid(-mid_rad:mid_rad, -mid_rad:mid_rad)
 sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)

 sel_pix_horiz = which(sel_grid[,2] <  35 & sel_grid[,2] > -35)
 sel_pix_vert = which(sel_grid[,1] <  35 & sel_grid[,1] > -35)

 good_star_cut$imDat[sel_pix_vert] = original[sel_pix_vert]
 good_star_cut$imDat[sel_pix_horiz] = original[sel_pix_horiz]

 }
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0 | is.na(good_star_cut_mask))) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  
  
  png(paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}

mid_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

mid_rad = (mid_scan[1,"dim_1"][[1]] -1L) / 2
sel_grid = expand.grid(-mid_rad:mid_rad, -mid_rad:mid_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 200 & sel_grid[,3] < 220)

mid_list =  foreach(ID = good_match_mid$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_im = Rfits_read(paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Mid_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.5){
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

stack_mid = propaneStackFlatFunc(image_list = mid_list, imager_func = parmed)

stack_mid_sd = propaneStackFlatFunc(image_list =  mid_list, imager_func = parsd)

rm(mid_list)
stack_mid$stdev = stack_mid_sd$image

Rfits_write(stack_mid, filename = paste0('/Path2StarFits/Mid_Star_Stack_',region,'_',band,'.fits'))


                                       ###########INNER STACK###########

dummy = foreach(ID = good_match_inner$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_inner$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_inner[ID,ra], GAIA_inner[ID,dec], box = 501, type='coord']
  
  png(paste0('/Path2StarFits/Inner_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Inner_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
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
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Inner_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  
  
  png(paste0('/Path2StarFits/Inner_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  
  return(NULL)
}


inner_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Inner_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

inner_rad = (inner_scan[1,"dim_1"][[1]] -1L) / 2 
sel_grid = expand.grid(-inner_rad:inner_rad, -inner_rad:inner_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 100 & sel_grid[,3] < 140)

inner_list = foreach(ID = good_match_inner$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_im = Rfits_read(paste0('/Path2StarFits/Inner_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Inner_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  if(temp_mask[[1]]$keyvalues$MASK < 0.02){
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

stack_inner = propaneStackFlatFunc(image_list = inner_list, imager_func = parmed)

stack_inner_sd = propaneStackFlatFunc(image_list =  inner_list, imager_func = parsd)

rm(inner_list)

stack_inner$stdev = stack_inner_sd$image

Rfits_write(stack_inner, filename = paste0('/Path2StarFits/Inner_Star_Stack_',region,'_',band,'.fits'))


                                       ###########CORE STACK###########

dummy = foreach(ID = good_match_core$bestmatch[,1])%dopar%{
  good_star = Rfits_read_image(temp_scan[good_match_core$ID[ID,1],"full"], ext=2)
  good_star_cut = good_star[GAIA_core[ID,ra], GAIA_core[ID,dec], box = 201, type='coord']
  
  png(paste0('/Path2StarFits/Core_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Core_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  
  segim = profoundProFound(good_star_cut, box=50)$segim #Run profound and select segmentation maps
 
  good_star_cut_mask = segim > 0 | is.na(good_star_cut_mask) #Mask where we have positive segim or NA
  
  good_star_cut$imDat[good_star_cut_mask > 0] = NA
  
  good_star_cut$keyvalues$MASK = length(which(good_star_cut_mask > 0)) / prod(dim(good_star_cut_mask))
  good_star_cut$keynames = c(good_star_cut$keynames, 'MASK')
  good_star_cut$keycomments = c(good_star_cut$keycomments, 'Fraction of masked pixels')
  

  Rfits_write_image(good_star_cut, paste0('/Path2StarFits/Core_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
  
  
  png(paste0('/Path2StarFits/Core_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.png'), width=1e3, height=1e3, type='cairo')
  plot(good_star_cut, qdiff=T)
  legend('topleft', legend=ID)
  dev.off()
  
  return(NULL)
}


core_scan = Rfits_key_scan(dirlist=paste0('/Path2StarFits/Core_Stars_',region,'_',band,'/'), get_dim=TRUE, cores=8)

core_rad = (core_scan[1,"dim_1"][[1]] -1L) / 2 
sel_grid = expand.grid(-core_rad:core_rad, -core_rad:core_rad)
sel_grid[,3] = sqrt(sel_grid[,1]^2 + sel_grid[,2]^2)
sel_pix = which(sel_grid[,3]> 0 & sel_grid[,3] < 100)


core_list = foreach(ID = good_match_core$bestmatch[,1], .errorhandling = 'remove')%dopar%{
  temp_im = Rfits_read(paste0('/Path2StarFits/Core_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_orig.fits'))
  temp_mask = Rfits_read(paste0('/Path2StarFits/Core_Stars_',region,'_',band,'/','GAIA_',formatC(ID,width=3,flag=0),'_mask.fits'))
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
    return(temp_im[[1]]$imDat)
  }else{
    return(stop('Skipping', ID))
  }
}

stack_core = propaneStackFlatFunc(image_list = core_list, imager_func = parmed)

stack_core_sd = propaneStackFlatFunc(image_list =  core_list, imager_func = parsd)

rm(core_list)

stack_core$stdev = stack_core_sd$image

Rfits_write(stack_core, filename = paste0('/Path2StarFits/Core_Star_Stack_',region,'_',band,'.fits'))
