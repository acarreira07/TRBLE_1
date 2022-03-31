############################
#####Initiate R session#####
############################
#Clear workspace
rm(list=ls());

#Set strings to be characters, not factors
options(stringsAsFactors=F);

#Load packages
library(glmnet);
library(caret);
library(MASS);
library(ape);
library(entropy);
library(markovchain);
library(igraph)
library(sBIC)

#Source TREBLE functions (change path to point R to wherever this is)
source("/Users/arnaldo/Desktop/TREBLE/treble_functions_220321.R")

###############################################################
#####Load functions (highlight this whole section and run)#####
###############################################################


#Bearing
bearing = function(x1=10, y1 = 10, x2=3, y2=3){
  require(NISTunits)
  if((x1 == x2) & (y1 > y2)){
    return(360)
  }else if((y1 == y2) & (x1 < x2)){
    return(90)
  } else if((y1==y2 & x1 > x2)){
    return(270)
  }else if(y1 == y2 & x1 < x2){
    return(180)
  }else if((x1 == x2) & (y1==y2)){
    return(NaN)
  }
  else
    theta = atan2(x2 - x1, y1 - y2)
  if(theta < 0){
    theta = theta + 2*pi
  }
  theta = NISTradianTOdeg(theta)
  return(theta)
}

#Function to calculate within condition variance
within_species_variance_umap = function(layout,
                                        extract_condition = FALSE,
                                        condition){
  
  if(extract_condition == TRUE){
    #Extract species
    l = layout[layout$strain%in%condition,]
  }else{
    l = layout
  }
  
  #Split on trial
  trials = split(l, l$trial)
  
  #Extract all possible bins
  xy_new_u = unique(layout$xy_new)
  
  #Create dataframe for results
  n = names(trials)
  
  bin_variance = as.data.frame(matrix(nrow = length(xy_new_u),
                                      ncol = length(n)))
  colnames(bin_variance) = n
  rownames(bin_variance) = xy_new_u
  
  for(j in 1:length(trials)){
    #print(j)
    t = table(trials[[j]]$xy_new)
    
    bin_variance[,j] = t[match(rownames(bin_variance),
                               names(t))]
  }
  
  bin_variance[is.na(bin_variance)] = 0
  
  #Calculate variance
  v = apply(bin_variance, 1, function(x) sd(x))
  
  s = colSums(bin_variance)
  
  bin_variance_percent = bin_variance
  for(i in 1:ncol(bin_variance)){
    bin_variance_percent[,i] = bin_variance_percent[,i]/s[i]
  }
  
  v_p = apply(bin_variance_percent, 1, function(x) sd(x))
  
  z = list(bin_variance,
           v,
           bin_variance_percent,
           v_p)
  names(z) = c("bin_counts_raw",
               "bin_variance_raw",
               "bin_counts_percent",
               "bin_variance_percent")
  return(z)
}

#Function to run glm using k-fold cross validation
run_glmnet_k_fold = function(outcome,
                             predictor_matrix,
                             family = 'gaussian',
                             train_prop = 0.8,
                             block_train = FALSE,
                             k = 10,
                             plot = FALSE){
  
  #Split data into test and train
  predictor_matrix$outcome = outcome
  if(block_train == TRUE){
    
    start = sample(seq(1, nrow(predictor_matrix)-round(nrow(predictor_matrix)*train_prop), 1), 1)
    train_inds = seq(start, start+round(nrow(predictor_matrix)*train_prop), 1)
    
  }else{
    train_inds = createDataPartition(y = predictor_matrix$outcome, 
                                     p = train_prop, 
                                     list = FALSE)
  }
  X_train = predictor_matrix[train_inds,]
  X_test = predictor_matrix[-train_inds,]
  
  #Fit model using 10-fold cross validation to get optimal alpha and lambda
  ctrl = trainControl(method = "cv", number = k)
  fit  =  train(outcome ~ ., 
                data = X_train,
                family = 'binomial', 
                method = "glmnet",
                preProc = c("center", "scale"),
                trControl = ctrl)
  #print(fit)
  
  #Get prediction
  yhats <- predict(fit, newdata = X_test, type = "raw")
  corr = cor(yhats, X_test$outcome)
  print(corr)
  
  if(plot == TRUE){
    plot(X_test$outcome,
         type = 'l',
         col = 'grey50',
         xlab = "Time",
         ylab = "Activity",
         cex.lab = 1.5,
         cex.axis = 1.5)
    lines(yhats,
          col = 'red')
  }
  
  #Return
  l = list(fit, yhats, X_train, X_test, corr)
  names(l) = c('mod', 'predicted_values', 'train_set', 'test_set', 'fit')
  return(l)
}

#Function to load full dataset
load_full_larvae_data = function(path_to_file){
  
  #Load
  x = read.csv(path_to_file, row.names = 1)
  colnames(x) = paste("larva_", seq(1, ncol(x), 1), sep = "")
  
  #Extract features
  x$feature = gsub("\\s*\\([^\\)]+\\)", "", rownames(x))
  x$time = gsub("[\\(\\)]", "", regmatches(rownames(x), gregexpr("\\(.*?\\)", rownames(x))))
  
  #Split on feature
  f = split(x, x$feature)
  
  #Combine into individual larvae
  trials = list()
  for(i in 1:(ncol(x)-2)){
    y = lapply(f, function(x) x[,i])
    
    y = as.data.frame(do.call(cbind, y))
    y$mom_theta = c(0, atan2(diff(as.numeric(y$mom_x)), diff(as.numeric(y$mom_y)))*(180/pi))
    y$head_theta = c(0, atan2(diff(as.numeric(y$head_x)), diff(as.numeric(y$head_y)))*(180/pi))
    y$tail_theta = c(0, atan2(diff(as.numeric(y$tail_x)), diff(as.numeric(y$tail_y)))*(180/pi))
    
    y$mom_vr = c(0, diff(y$mom_theta))
    y$head_vr = c(0, diff(y$head_theta))
    y$tail_vr = c(0, diff(y$tail_theta))
    
    y$time = seq(1, nrow(y), 1)
    
    trials[[colnames(x)[i]]] = y
  }
  
  return(trials)
}

#Function to calculate per bin entropies
calculate_bin_entropies = function(bin,
                                   strain_coords,
                                   strain_xy,
                                   window = 30,
                                   n_bins,
                                   plot_trajectories = FALSE,
                                   calculate_entropies = FALSE,
                                   n_trajectories = NULL){
  
  #Choose bin
  z = bin
  tmp = strain_coords
  
  #Split on bin and extract bouts in between
  y = which(!tmp == z)
  idx <- c(0, cumsum(abs(diff(y)) > 1))
  indices = split(y, idx)
  
  print("Extracting inter-bin bouts")
  #Get bin name
  series = lapply(indices, function(x) tmp[x])
  
  #Remove first bout (doesn't originate at bin of interest)
  series = series[-1]
  
  #Require bout of length n
  series = series[lapply(series, length)>=window]
  series = lapply(series, function(x) x[1:window])
  
  #Make empty vector for saving output
  entropies = c()
  
  if(length(series)>1){
    
    #Combine into dataframe
    series = do.call(cbind, series)
    print(ncol(series))
    
    #Add row corresponding to starting bin
    series = rbind(rep(z, ncol(series)), series)
    
    #Plot if desired
    if(plot_trajectories == TRUE){
      plot(series[,1],
           type = "l",
           lwd = 2,
           col = alpha("grey40", 0.5),
           ylim = c(min(unlist(series)), max(unlist(series))),
           ylab = "Bin",
           xlab = "Time",
           cex.axis = 1.5,
           cex.lab = 1.5,
           bty = 'n')
      for(i in 2:ncol(series)){
        lines(series[,i],
              lwd = 2,
              col = alpha("grey40", 0.5))}}
    
    if(!is.null(n_trajectories)){
      
      if(ncol(series)>n_trajectories){
        print("Calculating entropies")
        
        series = series[,sample(seq(1, ncol(series), 1), n_trajectories)]
        print(ncol(series))
        
        ##Calculate probability density functions
        for(h in 1:nrow(series)){
          
          #print(i)
          #Extract row(s) of interest
          row = unlist(series[h,])
          
          #Convert to coords
          row = cbind(as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[1]}))),
                      as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[2]}))))
          
          w = MASS::kde2d(row[,1], 
                          row[,2],
                          h = c(1, 1),
                          n = n_bins,
                          lims = c(1, n_bins+1, 1, n_bins+1))$z
          
          pdf = as.numeric(unlist(as.data.frame(w)))
          entropies = c(entropies, entropy::entropy(pdf, unit = 'log2'))
        }
        #plot(entropies,
        #type = "l",
        #lwd = 1.5,
        #cex.lab = 1.5,
        #cex.axis = 1.5,
        #xlab = "Time",
        #bty = 'n',
        #ylim = c(0, 4))
      }
      #else{
      #entropies = c(entropies, NA)
      #}
    }
    
    if(calculate_entropies == TRUE){
      print("Calculating entropies")
      ##Calculate probability density functions
      for(h in 1:nrow(series)){
        #print(i)
        #Extract row(s) of interest
        row = unlist(series[h,])
        
        #Convert to coords
        row = cbind(as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[1]}))),
                    as.numeric(unlist(lapply(strsplit(apply(as.data.frame(row), 1, function(x) unique(strain_xy[strain_coords == x])), "_"), function(v){v[2]}))))
        
        w = MASS::kde2d(row[,1], 
                        row[,2],
                        h = c(1, 1),
                        n = n_bins,
                        lims = c(1, n_bins+1, 1, n_bins+1))$z
        
        pdf = as.numeric(unlist(as.data.frame(w)))
        entropies = c(entropies, entropy::entropy(pdf, unit = 'log2'))
      }
      #plot(entropies,
      #type = "l",
      #lwd = 1.5,
      #cex.lab = 1.5,
      #cex.axis = 1.5,
      #xlab = "Time",
      #bty = 'n',
      #ylim = c(0, 4))
    }
    
    l = list(series, entropies)
    names(l) = c("bouts", "entropies")
    return(l)
  }
}

#Markov generative test
markov_generative_test = function(trials,
                                  states,
                                  reps = 10,
                                  n_sequences = 1000,
                                  sequence_length = 5,
                                  other_to_compare = NULL,
                                  random_transitions = FALSE){
  MC = markovchain::markovchainFit(data=states)$estimate
  
  #Replace transition matrix with random values if desired
  if(random_transitions == TRUE){
    r = matrix(ncol = ncol(MC@transitionMatrix), nrow = 0)
    for(i in 1:nrow(MC@transitionMatrix)){
      
      t = runif(ncol(MC@transitionMatrix))
      t = t/sum(t)
      r = rbind(r, t)
    }
    
    colnames(r) = colnames(MC@transitionMatrix)
    rownames(r) = rownames(MC@transitionMatrix)
    
    MC@transitionMatrix = r
  }
  
  print('run generative test')
  #Generative test
  f = c()
  frs = c()
  fss = c()
  
  for(h in 1:reps){
    
    Sys.sleep(.1)
    cat(h)
    flush.console()
    
    #Generate n number of synthetic sequences
    n = n_sequences
    x = sequence_length
    seqs = c()
    for(i in 1:n){
      seqs = c(seqs, paste(markovchainSequence(n=x, markovchain=MC), collapse = ""))
    }
    
    #Create library of sequences to match to
    #Randomly sample n sequences of length x from each trial structured by trial length
    lengths = sapply(trials, function(x) length(as.character(na.omit(x$louvain_cluster))))
    lengths = lengths/sum(lengths)
    
    exp = round(lengths*n)
    rest = n-sum(exp)
    
    #If expectation is too large or lower adjust
    if(rest>0){
      exp[order(exp, decreasing = TRUE)[1]] = exp[order(exp, decreasing = TRUE)[1]]+rest
    }
    if(rest<0){
      exp[order(exp, decreasing = TRUE)[1]] = exp[order(exp, decreasing = TRUE)[1]]+rest
    }
    
    #Create library dfs
    lib1 = c()
    if(is.null(other_to_compare) == FALSE){
      lengths = sapply(other_to_compare, function(x) length(as.character(na.omit(x$louvain_cluster))))
      lengths = lengths/sum(lengths)
      
      exp2 = round(lengths*n)
      
      for(i in 1:length(other_to_compare)){
        
        #print(i)
        #Get vector of just behaviors
        y = as.character(na.omit(other_to_compare[[i]]$louvain_cluster))
        
        #Randomly sample n seqs of length x
        ints = sample(seq(1, length(y)-x, 1), exp2[i])
        
        #Grab sequences
        for(j in 1:length(ints)){
          z = paste(y[ints[j]:(ints[j]+(x-1))], collapse = "")
          lib1 = rbind(lib1, z)
        }
      }
    }else{
      for(i in 1:length(trials)){
        
        #print(i)
        #Get vector of just behaviors
        y = as.character(na.omit(trials[[i]]$louvain_cluster))
        
        #Randomly sample n seqs of length x
        ints = sample(seq(1, length(y)-x, 1), exp[i])
        
        #Grab sequences
        for(j in 1:length(ints)){
          z = paste(y[ints[j]:(ints[j]+(x-1))], collapse = "")
          lib1 = rbind(lib1, z)
        }
      }
    }
    
    lib2 = c()
    for(i in 1:length(trials)){
      
      #Get vector of just behaviors
      y = as.character(na.omit(trials[[i]]$louvain_cluster))
      
      #Randomly sample n seqs of length x
      ints = sample(seq(1, length(y)-x, 1), exp[i])
      
      #Grab sequences
      for(j in 1:length(ints)){
        z = paste(y[ints[j]:(ints[j]+(x-1))], collapse = "")
        lib2 = rbind(lib2, z)
      }
    }
    
    #Compare synthetic to real
    fs = table(lib1%in%seqs)
    
    if(is.na(fs[2])){
      fs[2] = 0
    }
    
    fss = c(fss, as.numeric(fs[2]))
    
    #Compare real to real
    fr = table(lib1%in%lib2)
    
    if(is.na(fr[2])){
      fr[2] = 0
    }
    
    frs = c(frs, as.numeric(fr[2]))
    
    #Calculate f
    f = c(f, as.numeric(fs[2])/as.numeric(fr[2]))
  }
  
  l = list(fss, frs, f, MC)
  names(l) = c("f_synthetic", "f_real", "f", "markov_model")
  return(l)
}


##################################################
#####Load and clean data (collected at 10hz)######
##################################################
##In this analysis we were comparing larvae under light and no light conditions 
#####Loading no light#####
#Initiate empty list to load files into
no_light = list()

###&&&Change_path&&&###

#Set working directory to where the behavior files are (you'll need to change the path/name to match your data)
setwd('/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_no_light')

#List files in directory
files = list.files()

#Loop through and load files into 'no_light' list
for(i in 1:length(files)){
  print(paste(i, 'out of', length(files)))
  no_light[[i]] = load_full_larvae_data(files[i])
}

#Combine
no_light = do.call(c, no_light)

#Change names and add larvae
names(no_light) = paste(rep('no_light', length(no_light)), seq(1, length(no_light), 1), sep = '_')
for(i in 1:length(no_light)){
  
  #Remove NA rows
  no_light[[i]] = no_light[[i]][rowSums(is.na(no_light[[i]])) != ncol(no_light[[i]]),]
  
  #Remove NAs
  no_light[[i]] = no_light[[i]][complete.cases(no_light[[i]]),]
  
  #Add larvae
  no_light[[i]]$larvae = rep(names(no_light)[i], nrow(no_light[[i]]))
}

#Turn into matrix
no_light = do.call(rbind, no_light)

#####Loading light#####
#Initiate empty list to load files into
light = list()


###&&&Change_path&&&###


#Set working directory to where the behavior files are (you'll need to change the path/name to match your data)
setwd('/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_plus_light')

#List files in directory
files = list.files()

#Loop through and load files into 'light' list
for(i in 1:length(files)){
  print(paste(i, 'out of', length(files)))
  light[[i]] = load_full_larvae_data(files[i])
}

#####Combining light and no_light into one matrix and cleaning#####
#Combine
light = do.call(c, light)

#Change names and add larvae
names(light) = paste(rep('light', length(light)), seq(1, length(light), 1), sep = '_')
for(i in 1:length(light)){
  
  #Remove NA rows
  light[[i]] = light[[i]][rowSums(is.na(light[[i]])) != ncol(light[[i]]),]
  
  #Remove NAs
  light[[i]] = light[[i]][complete.cases(light[[i]]),]
  
  #Add larvae
  light[[i]]$larvae = rep(names(light)[i], nrow(light[[i]]))
}

#Turn into matrix
light = do.call(rbind, light)

#Combine light and no_light into one matrix
y3 = as.data.frame(do.call(rbind, list(no_light, light)))

#Remove NA rows
y3 = y3[rowSums(is.na(y3)) != ncol(y3),]

#Remove NAs
y3 = y3[complete.cases(y3),]

#Remove not well oriented rows?
y3 = y3[!y3$is_well_oriented == 0,]
y3 = y3[!y3$is_coiled == 1,]

#Convert to proper data types
y3 = as.data.frame(y3)
for(i in 1:37){
  y3[,i] = as.numeric(y3[,i])
}

#Filter to desired features
y = data.frame(area = y3$area,
               #bending = y3$bending,
               bending = abs((y3$bending-180)),
               velocity = y3$velocity,
               spine = y3$spine_length,
               radius_1 = y3$radius_1,
               radius_2 = y3$radius_2,
               radius_3 = y3$radius_3,
               perimeter = y3$perimeter,
               head_vr = abs(y3$head_vr),
               mom_vr = abs(y3$mom_vr),
               tail_vr = abs(y3$tail_vr),
               dist = y3$dst_to_origin,
               row.names = rownames(y3))

#Remove NA rows
y = y[rowSums(is.na(y)) != ncol(y),]

#Remove NAs
y = y[complete.cases(y),]

#Split
s = split(y, 
          unlist(lapply(strsplit(rownames(y), "\\."), function(v){v[1]})))

###&&&Change_path&&&###

#Save (change the path, etc. to match your needs)
saveRDS(s, '/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_features_individual_trials_220329.RDS')

#####Filtering and detrending, if desired#####
##Here, we required trials to have at  least 250 timepoints
#Select only trials longer than n (if desired)
s = s[lapply(s, function(x) nrow(x))>250]

#Filter on distance traveled (if desired)
s = s[lapply(s, function(x) max(x$dist, na.rm = TRUE))>50]

#Detrend size values (if desired; generally recommended since, as larvae get further from center, their recorded size varies)
for(i in 1:length(s)){
  print(i)
  s[[i]]$area = s[[i]]$area/forecast::ma(s[[i]]$area,10)
  s[[i]]$area[is.na(s[[i]]$area)] = median(s[[i]]$area, na.rm = TRUE)
  s[[i]]$perimeter = s[[i]]$perimeter/forecast::ma(s[[i]]$perimeter,10)
  s[[i]]$perimeter[is.na(s[[i]]$perimeter)] = median(s[[i]]$perimeter, na.rm = TRUE)
  s[[i]]$radius_1 = s[[i]]$radius_1/forecast::ma(s[[i]]$radius_1,10)
  s[[i]]$radius_1[is.na(s[[i]]$radius_1)] = median(s[[i]]$radius_1, na.rm = TRUE)
  s[[i]]$radius_2 = s[[i]]$radius_2/forecast::ma(s[[i]]$radius_2,10)
  s[[i]]$radius_2[is.na(s[[i]]$radius_2)] = median(s[[i]]$radius_2, na.rm = TRUE)
  s[[i]]$radius_3 = s[[i]]$radius_3/forecast::ma(s[[i]]$radius_3,10)
  s[[i]]$radius_3[is.na(s[[i]]$radius_3)] = median(s[[i]]$radius_3, na.rm = TRUE)
  s[[i]]$spine = s[[i]]$spine/forecast::ma(s[[i]]$spine,10)
  s[[i]]$spine[is.na(s[[i]]$spine)] = median(s[[i]]$spine, na.rm = TRUE)
}

#Calculate z-scores
for(i in 1:length(s)){
  for(j in 1:11){
    s[[i]][,j] = scale(s[[i]][,j])
    #s[[i]][,j] = s[[i]][,j]/max(s[[i]][,j])
  }
}

#Recombine
y = do.call(rbind, s)

###&&&Change_path&&&###

#Save (change the path, etc. to match your needs)
saveRDS(y, '/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_features_210329.RDS')

#####Running PCA on input  parameters#####
##Since some of the behavioral parameters that were tracked are correlated, PCA can help reduce the number of dimensions needed
#Run PCA
pca = prcomp(y[,!colnames(y) == 'dist'], scale. = TRUE, center = TRUE)
summary(pca)

#Plot
plot(as.matrix(summary(pca)[[6]])[3,]*100, type = 'b', 
     pch = 20,
     xlab = 'n PCs',
     ylab = 'Variation explained',
     ylim = c(0,100),
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex = 2,
     bty = 'n')
abline(h = 90,
       lty = 'dashed',
       col = 'gray60')

#Get pca loadings
##The first 8 pcs tend to capture the majority of variation in larval input parameters
s = pca$x[,1:8]

###&&&Change_path&&&###

#Save (change the path, etc. to match your needs)
saveRDS(s, '/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_PCA_220329.RDS')

#Clean
rm(y1, y2, y3, y, light, no_light, i, j, files, res, pca)
gc()

##########################################################
#####Get feature windows of desired size and run UMAP#####
##########################################################

# this part takes about 40min to run
#Get windows (here using a window size of 8, recommended for larval behavior)
win = get_windows(s, window_size = 8)

#Run umap
u = umap(t(win), verbose = TRUE)

#Extract UMAP 2d layout
layout = data.frame(x = u$layout[,1],
                    y = u$layout[,2])

#Bin (64x64 grid ('n_bins' option) recommended)
layout = bin_umap(layout, n_bins = 64)$layout

#Add trial
layout$trial = unlist(lapply(strsplit(rownames(layout), "\\."), function(v){v[1]}))

#Add time
layout$time = as.numeric(unlist(lapply(strsplit(rownames(layout), "\\."), function(v){v[3]})))

#Add id
layout$id = paste(layout$trial, layout$time, sep = '_')

#Add condition
layout$condition = unlist(lapply(strsplit(layout$trial, "_"), function(v){v[1]}))

###&&&Change_path&&&###

#Save (change to match your needs)
saveRDS(layout,
        '/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_windowsize_8_220329.RDS')

#######################################
#####Analyzing UMAP behavior space#####
#######################################

## For some reason i had to run this twice. 
## First time it did not plot the density plot and  gave an error. 
###&&&Change_path&&&###

#Load layout (wherever and with whichever name you decide to use)
layout = readRDS('/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_windowsize_8_220329.RDS')

#Plot 2d structure of UMAP behavior space
plot(layout[,1:2], 
     type = 'l',
     lwd = 0.1, 
     col = alpha('grey50', 0.1),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '',
     bty = 'n')

#Plot 2d probability density function of behavior space
image(kde2d(layout$x, 
            layout$y, 
            h = 1, 
            n = 300),
      col = c('white', hcl.colors(12, "YlOrRd", rev = TRUE)),
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n',
      xlab = '',
      ylab = '')

#Plot as a vector field
plot_vector_field(layout, bin_umap = TRUE, n_bins = 16, color_by_theta = TRUE)

##############################################################################################
#####Analyzing input features (e.g. radius, spine length) as a function of behavior space#####
##############################################################################################


###&&&Change_path&&&###


#Load features (wherever and with whichever name you decide to use)
features = readRDS('/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_features_210329.RDS')

#Match features to the rownames of the layout file (i.e. making sure they represent the same timepoints)
features = features[match(rownames(layout), rownames(features)),]

#Combine
layout = cbind(layout, features)

#Correlation between layout xy position and features
cor(layout$x, layout[,11:22])
cor(layout$y, layout[,11:22])

#Ploting features as means of layout bins
par(mfrow = c(2,6), mar = c(1,1,1,1))
x = split(layout, layout$xy_new)
for(i in 11:22){
  #z = round(layout[,i], 2)
  z = round(unlist(lapply(x, function(y) mean(y[[i]]))), 2)
  z[z>quantile(z, probs = 0.99)] = quantile(z, probs = 0.99)
  z[z<quantile(z, probs = 0.01)] = quantile(z, probs = 0.01)
  
  if(min(z, na.rm = TRUE)>0|max(z, na.rm = TRUE)<0){
    
    cols = colorRampPalette(c('grey90', 'darkred'))(length(seq(min(z, na.rm = TRUE), max(z, na.rm = TRUE), 0.01)))
    names(cols) = round(seq(min(z), max(z), 0.01), 2)
    cols = cols[match(z, names(cols))]
    
  }else{
    
    cols1 = colorRampPalette(c('midnightblue', 'grey90'))(length(seq(min(z, na.rm = TRUE), 0, 0.01)))
    cols2 = colorRampPalette(c('grey90', 'darkred'))(length(seq(0.01, max(z, na.rm = TRUE), 0.01)))
    
    names(cols1) = round(seq(min(z), 0, 0.01), 2)
    names(cols2) = round(seq(0.01, max(z), 0.01), 2)
    
    cols = c(cols1, cols2)
    cols = cols[match(z, names(cols))]
  }
  
  # plot(layout$x,
  #      layout$y,
  #      pch = 20,
  #      cex = 0.2,
  #      col = cols,
  #      bty = 'n',
  #      xaxt = 'n',
  #      yaxt = 'n',
  #      xlab = '',
  #      ylab = '')
  
  plot(unlist(lapply(strsplit(names(z), "_"), function(v){v[1]})),
       unlist(lapply(strsplit(names(z), "_"), function(v){v[2]})),
       pch = 20,
       #cex = 0.75,
       col = cols,
       bty = 'n',
       xaxt = 'n',
       yaxt = 'n',
       xlab = '',
       ylab = '')
  
  title(main = colnames(layout)[i],
        cex.main = 1.5,
        font.main = 1)
  
}

###&&&Change_path&&&###

#Save annotated (change to match your needs)
saveRDS(layout, '/Users/arnaldo/Desktop/TREBLE/mz1407/mz1407_umap_layout_220329.RDS')

###############################################################################
#####Comparing light and no_light conditions via density in behavior space#####
###############################################################################
#Split layout based on condition (to get a list with 2 elements: light layout and no_light layout)
condition = split(layout, layout$condition)

#Calculate the probability density function for the light condition
light_pdf = kde2d(condition$light$x, condition$light$y, h = 1, n = 300,
                  lims = c(c(min(layout$x), max(layout$x)),
                           c(min(layout$y), max(layout$y))))$z
light_pdf = light_pdf/max(light_pdf)

#Calculate the probability density function for the no_light condition
no_pdf = kde2d(condition$no$x, condition$no$y, h = 1, n = 300,
               lims = c(c(min(layout$x), max(layout$x)),
                        c(min(layout$y), max(layout$y))))$z
no_pdf = no_pdf/max(no_pdf)

#Calculate difference between the two probability density functions
d = light_pdf-no_pdf
d[d<(-0.15)] = -0.15
d[d>0.15] = 0.15

#Plot
par(mfrow = c(1,3))
image(light_pdf, bty = 'n', xaxt = 'n', yaxt = 'n', col = c('white', hcl.colors(12, "Greens", rev = TRUE)))
image(no_pdf, bty = 'n', xaxt = 'n', yaxt = 'n', col = c('white', hcl.colors(12, "BuPu", rev = TRUE)))
image(d, 
      col = colorRampPalette(c('darkmagenta', 'white', 'darkgreen'))(200), 
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n')
