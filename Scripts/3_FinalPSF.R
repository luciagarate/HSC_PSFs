#We first load the required R packages:
library(Rfits)
library(doParallel)
library(foreach)
library(ProPane)
library(imager) #parmed function in propaneStackFlatFunc

registerDoParallel(cores=8)

#Per band, we stack the four PSFs that we have (one per each GAMA-field) to create the final HSC-PDR3 PSF:
band = 'g'
region = c('G02','G09','G12','G15')

#Per band, we create a list with the four PSFs (the PSFs from each GAMA-field are already normalised):
stack_list = foreach(i = region, .errorhandling = 'remove')%dopar%{
  image = Rfits_read(paste0('/Path2PSFs/Combined_Star_Stack_',band,'_',i,'.fits')) 
  return(image[[1]]$imDat)
}

#We stack them:
stack_final = propaneStackFlatFunc(image_list = stack_list, imager_func = parmed)

stack_final_sd = propaneStackFlatFunc(image_list =  stack_list, imager_func = parsd)

stack_final$stdev = stack_final_sd$image

#We can visualise the PSF:
library(magicaxis)
magimage(stack_final$image, qdiff=T)

#Final PSF:
Rfits_write(stack_final, filename = paste0('/Path2FinalPSFs/Final_Star_Stack_',band,'.fits')) #Add compress = T if we want a compressed version