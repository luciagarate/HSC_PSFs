library(data.table)
library(celestial)
library(Rfits)
library(ProFound)
library(doParallel)
library(foreach)
#library(imager)


region = 'G02'
band = 'g'

#In this script we combine the outer, middle, inner and core stacks to create a final PSF per GAMA region and band
stack_bright = Rfits_read(paste0('/Path2StarFits/Bright_Star_Stack_',region,'_',band,'.fits'), pointer=FALSE) #outer stack
stack_combine = stack_bright$image$imDat

#We expand the grid with x,y positions of each pixel by adding a third column with the radius
bright_rad = (dim(stack_bright$image)[1] - 1L)/2 
sel_grid_b = expand.grid(-bright_rad:bright_rad, -bright_rad:bright_rad)
sel_grid_b[,3] = sqrt(sel_grid_b[,1]^2 + sel_grid_b[,2]^2)
sel_pix_b = which(sel_grid_b[,3]> 150 & sel_grid_b[,3] < 160) #this is the normalisation annulus that we use to combine the
#outer stack and the middle stack
scale_b = median(stack_combine[sel_pix_b], na.rm = TRUE)
stack_combine = stack_combine/scale_b #normalise the outer stack

stack_mod = Rfits_read(paste0('/Path2StarFits/Mod_Star_Stack_',region,'_',band,'.fits'), pointer=FALSE) #middle stack

mod_rad = (dim(stack_mod$image)[1] - 1L)/2 
sel_grid_m = expand.grid(-mod_rad:mod_rad, -mod_rad:mod_rad)
sel_grid_m[,3] = sqrt(sel_grid_m[,1]^2 + sel_grid_m[,2]^2)
sel_pix_m = which(sel_grid_m[,3]> 150 & sel_grid_m[,3] < 160)
scale_m = median(stack_mod$image$imDat[sel_pix_m], na.rm = TRUE)
stack_mod$image$imDat = stack_mod$image$imDat/scale_m #normalise the middle stack


#We now replace the pixels with r<150 in the outer stack with the middle stack:
fill_b = which(sel_grid_b[,3] < 150)
fill_m = which(sel_grid_m[,3] < 150)

stack_combine[fill_b] = stack_mod$image$imDat[fill_m]  #outer + middle combination


#We now combine the outer+ middle combination with the inner stack:
stack_faint = Rfits_read(paste0('/Path2StarFits/Faint_Star_Stack_',region,'_',band,'.fits'), pointer=FALSE) #inner stack

bright_rad = (dim(stack_bright$image)[1] - 1L)/2
sel_grid_b = expand.grid(-bright_rad:bright_rad, -bright_rad:bright_rad)
sel_grid_b[,3] = sqrt(sel_grid_b[,1]^2 + sel_grid_b[,2]^2)
sel_pix_b = which(sel_grid_b[,3]> 60 & sel_grid_b[,3] < 70) #we normalise with this annulus the outer+ middle and the inner
scale_b = median(stack_combine[sel_pix_b], na.rm = TRUE)
stack_combine = stack_combine/scale_b

faint_rad = (dim(stack_faint$image)[1] - 1L)/2
sel_grid_f = expand.grid(-faint_rad:faint_rad, -faint_rad:faint_rad)
sel_grid_f[,3] = sqrt(sel_grid_f[,1]^2 + sel_grid_f[,2]^2)
sel_pix_f = which(sel_grid_f[,3]> 60 & sel_grid_f[,3] < 70)
scale_f = median(stack_faint$image$imDat[sel_pix_f], na.rm = TRUE)
stack_faint$image$imDat = stack_faint$image$imDat/scale_f

fill_b = which(sel_grid_b[,3] < 60)
fill_f = which(sel_grid_f[,3] < 60)

#We replace the pixels with r<60 in the outer+middle combination with the inner stack:
stack_combine[fill_b] = stack_faint$image$imDat[fill_f] #outer+middle+inner combination

#We now combine the outer+ middle+inner combination with the core stack:
stack_faintest = Rfits_read(paste0('/Path2StarFits/Faintest_Star_Stack_',region,'_',band,'.fits'), pointer=FALSE)  #core stack

bright_rad = (dim(stack_bright$image)[1] - 1L)/2
sel_grid_b = expand.grid(-bright_rad:bright_rad, -bright_rad:bright_rad)
sel_grid_b[,3] = sqrt(sel_grid_b[,1]^2 + sel_grid_b[,2]^2)
sel_pix_b = which(sel_grid_b[,3]> 20 & sel_grid_b[,3] < 30) #normalise the outer+ middle+inner combination and the core stack w/this annulus
scale_b = median(stack_combine[sel_pix_b], na.rm = TRUE)
stack_combine = stack_combine/scale_b

faintest_rad = (dim(stack_faintest$image)[1] - 1L)/2
sel_grid_ff = expand.grid(-faintest_rad:faintest_rad, -faintest_rad:faintest_rad) #entre -100 y 100
sel_grid_ff[,3] = sqrt(sel_grid_ff[,1]^2 + sel_grid_ff[,2]^2)
sel_pix_ff = which(sel_grid_ff[,3]> 20 & sel_grid_ff[,3] < 30)
scale_ff = median(stack_faintest$image$imDat[sel_pix_ff], na.rm = TRUE)
stack_faintest$image$imDat = stack_faintest$image$imDat/scale_ff

#Replace the pixels with r<20 in the outer+middle+inner combination with the core stack:
fill_b = which(sel_grid_b[,3] < 20)
fill_ff = which(sel_grid_ff[,3] < 20)

stack_combine[fill_b] = stack_faintest$image$imDat[fill_ff] #final combination = PSF

#A few final details:
#Apply a symmetric process by making horizontal, vertical an diagonal reflections to obtain a symmetric PSF:
stack_combine = stack_combine + stack_combine[,4001:1] + stack_combine[4001:1,] + stack_combine[4001:1,4001:1]

stack_combine/sum(stack_combine) #normalise to 1 the flux of each PSF

Rfits_write_image(stack_combine, filename = paste0('/Path2StarFits/Combined_Star_Stack_',band,'_',region,'.fits'))