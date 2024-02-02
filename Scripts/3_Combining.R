#We first load the required R packages:
library(Rfits)
library(doParallel)
library(foreach)
library(ProPane)
library(imager) #parmed function in propaneStackFlatFunc


#In this script we combine the outer, middle, inner and core stacks to create a final PSF per band
stack_outer = Rfits_read(paste0('/Path2StarFits/Outer_Star_Stack_',band,'.fits'), pointer=FALSE) #outer stack
stack_combine = stack_outer$image$imDat

#We expand the grid with x,y positions of each pixel by adding a third column with the radius
outer_rad = (dim(stack_outer$image)[1] - 1L)/2 
sel_grid_o = expand.grid(-outer_rad:outer_rad, -outer_rad:outer_rad)
sel_grid_o[,3] = sqrt(sel_grid_o[,1]^2 + sel_grid_o[,2]^2)
sel_pix_o = which(sel_grid_o[,3]> 150 & sel_grid_o[,3] < 160) #this is the normalisation annulus that we use to combine the outer stack and the middle stack
scale_o = median(stack_combine[sel_pix_o], na.rm = TRUE)
stack_combine = stack_combine/scale_o #normalise the outer stack

stack_mid = Rfits_read(paste0('/Path2StarFits/Mid_Star_Stack_',band,'.fits'), pointer=FALSE) #middle stack

mid_rad = (dim(stack_mid$image)[1] - 1L)/2 
sel_grid_m = expand.grid(-mid_rad:mid_rad, -mid_rad:mid_rad)
sel_grid_m[,3] = sqrt(sel_grid_m[,1]^2 + sel_grid_m[,2]^2)
sel_pix_m = which(sel_grid_m[,3]> 150 & sel_grid_m[,3] < 160)
scale_m = median(stack_mid$image$imDat[sel_pix_m], na.rm = TRUE)
stack_mid$image$imDat = stack_mid$image$imDat/scale_m #normalise the middle stack


#We now replace the pixels with R<150 in the outer stack with the middle stack:
fill_o = which(sel_grid_o[,3] < 150)
fill_m = which(sel_grid_m[,3] < 150)

stack_combine[fill_o] = stack_mid$image$imDat[fill_m]  #outer + middle combination


#We now combine the outer + middle combination with the inner stack:
stack_inner = Rfits_read(paste0('/Path2StarFits/Inner_Star_Stack_',band,'.fits'), pointer=FALSE) #inner stack

outer_rad = (dim(stack_outer$image)[1] - 1L)/2
sel_grid_o = expand.grid(-outer_rad:outer_rad, -outer_rad:outer_rad)
sel_grid_o[,3] = sqrt(sel_grid_o[,1]^2 + sel_grid_o[,2]^2)
sel_pix_o = which(sel_grid_o[,3]> 60 & sel_grid_o[,3] < 70) #we normalise with this annulus the outer+ middle and the inner
scale_o = median(stack_combine[sel_pix_o], na.rm = TRUE)
stack_combine = stack_combine/scale_o

inner_rad = (dim(stack_inner$image)[1] - 1L)/2
sel_grid_i = expand.grid(-inner_rad:inner_rad, -inner_rad:inner_rad)
sel_grid_i[,3] = sqrt(sel_grid_i[,1]^2 + sel_grid_i[,2]^2)
sel_pix_i = which(sel_grid_i[,3]> 60 & sel_grid_i[,3] < 70)
scale_i = median(stack_inner$image$imDat[sel_pix_i], na.rm = TRUE)
stack_inner$image$imDat = stack_inner$image$imDat/scale_i

fill_o = which(sel_grid_o[,3] < 60)
fill_i = which(sel_grid_i[,3] < 60)

#We replace the pixels with r<60 in the outer+middle combination with the inner stack:
stack_combine[fill_o] = stack_inner$image$imDat[fill_i] #outer+middle+inner combination

#We now combine the outer+ middle+inner combination with the core stack:
stack_core = Rfits_read(paste0('/Path2StarFits/Core_Star_Stack_',band,'.fits'), pointer=FALSE)  #core stack

outer_rad = (dim(stack_outer$image)[1] - 1L)/2
sel_grid_o = expand.grid(-outer_rad:outer_rad, -outer_rad:outer_rad)
sel_grid_o[,3] = sqrt(sel_grid_o[,1]^2 + sel_grid_o[,2]^2)
sel_pix_o = which(sel_grid_o[,3]> 5 & sel_grid_o[,3] < 30) #normalise the outer+ middle+inner combination and the core stack w/this annulus
scale_o = median(stack_combine[sel_pix_o], na.rm = TRUE)
stack_combine = stack_combine/scale_o

core_rad = (dim(stack_core$image)[1] - 1L)/2
sel_grid_c = expand.grid(-core_rad:core_rad, -core_rad:core_rad) #entre -100 y 100
sel_grid_c[,3] = sqrt(sel_grid_c[,1]^2 + sel_grid_c[,2]^2)
sel_pix_c = which(sel_grid_c[,3]> 5 & sel_grid_c[,3] < 30)
scale_c = median(stack_core$image$imDat[sel_pix_c], na.rm = TRUE)
stack_core$image$imDat = stack_core$image$imDat/scale_c

#Replace the pixels with r<20 in the outer+middle+inner combination with the core stack:
fill_o = which(sel_grid_o[,3] < 20)
fill_c = which(sel_grid_c[,3] < 20)

stack_combine[fill_o] = stack_core$image$imDat[fill_c] #final combination = PSF

#A few final details:
#Apply a symmetric process by making horizontal, vertical an diagonal reflections to obtain a symmetric PSF: (comment if the user does not want to symmetrise the PSF)
stack_combine = stack_combine + stack_combine[,4001:1] + stack_combine[4001:1,] + stack_combine[4001:1,4001:1]

stack_combine/sum(stack_combine) #normalise to 1 the flux of each PSF

Rfits_write_image(stack_combine, filename = paste0('/Path2PSFs/HSCPSF_',band,'.fits'))
