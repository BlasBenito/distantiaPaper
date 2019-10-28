#generates combinations of parameters to feed a sensitivity analysis
#data: dataframe, required to know the dimensions of the input data.
generateParams <- function(data, transformations){
  
  #generating combinations of parameters
  params <- expand.grid(
    method = c("lock-step", "elastic"), 
    transformation = transformations, 
    distance = c("euclidean", "manhattan"), 
    ignore.blocks = c(TRUE, FALSE),
    diagonal = c(TRUE, FALSE), 
    paired.samples = c(TRUE, FALSE),
    cases.to.modify = c(10, 100, 200, 400, 800, 1600, 3200), 
    remove.samples = floor(
      seq(
      0, 
      floor(nrow(data) - (nrow(data) / 10)), 
      length.out = 7
    )
    )
    , 
    stringsAsFactors = FALSE
    )
  
  #removing exceptions
  params$remove <- "NO"
  #lock-step doesn't have diagonal option
  #we remove half of them
  params[params$method == "lock-step" & params$diagonal == TRUE, "remove"] <- NA
  
  #sequences must have same number of cases
  params[params$method == "lock-step" & params$remove.samples != 0, "remove"] <- NA
  
  #lock-step requires paired samples
  params[params$method == "lock-step" & params$paired.samples == FALSE, "remove"] <- NA
  
  #lock-step does not require ignore.blocks
  params[params$method == "lock-step" & params$ignore.blocks == TRUE, "remove"] <- NA
  
  #elastic mode doesn't need paired samples
  params[params$method == "elastic" & params$paired.samples == TRUE, "remove"] <- NA
  
  #we are not exploring the interaction between cases.to.modify and remove.samples
  params <- rbind(
    params[params$remove.samples == 0,], 
    params[params$cases.to.modify == 0,]
  )
  
  #removing NA
  params <- na.omit(params)
  
  #removing duplicates
  params <- params[!duplicated(params), ]
  
  #removing a column
  params$remove <- NULL
  
  #ordering
  params <- params[with(params, order(method, transformation, distance, diagonal)), ]
  
  return(params)
  
}


#x: dataframe, sequence to analyze.
#transformation: char, name of the transformation.
#distance: char, name of the distance metric.
#diagonal: boolean, switches diagonal in elastic mode.
#paired.samples: boolean, TRUE for lock-step mode.
#cases.to.modify: integer, number of data-points to change.
#remove.samples: integer (< nrow(x)), cases to remove.
#ignore.blocks: boolean, ignores blocks if TRUE
#repetitions: integer, number of repetitions.
repetitionsPsi <- function(
  x,
  transformation,
  distance,
  diagonal,
  paired.samples,
  cases.to.modify,
  remove.samples,
  ignore.blocks,
  repetitions
){
  
  require(distantia)
  
  #getting only numeric values of x
  x <- x[, unlist(lapply(x, is.numeric))]
  
  #psi values
  psi.values <- rep(NA, repetitions)
  
  #iterating through repetitions
  for(i in 1:repetitions){
    
    #creating a copy
    x.copy <- x
    
    #removing samples
    if(remove.samples > 0){
      x.copy <- x.copy[-sample(nrow(x.copy), remove.samples), ]
    }
    
    
    #iterating through cases.to.modify
    if(cases.to.modify > 0){
      
      #select columns
      selected.cols <- sample(ncol(x), cases.to.modify, replace = TRUE)
      
      #select rows
      selected.rows <- sample(nrow(x), cases.to.modify, replace = TRUE)
      
      #changing cases
      for(j in 1:cases.to.modify){
        
        val <- x.copy[selected.rows[j], selected.cols[j]]
        
        #changing value of the cell by adding or substracting a fraction of it (between 1 and 100%)
        x.copy[selected.rows[j], selected.cols[j]] <- val + sample(c(1, -1), 1) * (val / sample(1:100, 1))
        
        #changing case by adding a 10% of its value
        # x.copy[selected.rows[j], selected.cols[j]] <- val + (val / 10)
      }
      
    }
    
    #preparing sequences
    #transformation = "none"
    x.joint <- distantia::prepareSequences(
      sequence.A = x,
      sequence.A.name = "x",
      sequence.B = x.copy,
      sequence.B.name = "x.copy",
      grouping.column = "id",
      transformation = transformation
    )
    
    #psi transformation.1 - manhattan
    psi.values[i] <- workflowPsi(
      sequences = x.joint,
      grouping.column = "id",
      method = distance,
      diagonal = diagonal,
      paired.samples = paired.samples,
      ignore.blocks = ignore.blocks,
      parallel.execution = FALSE
    )$psi
    
  }
  
  return(psi.values)
  
}


#RUNS THE SENSITIVITY ANALYSIS
#data: dataframe with the target sequence
#params: dataframe with parameters
#repetitions: integer, number of repetitions
#to pass to repetitionsPsi
sensitivityPsi <- function(data, params, repetitions){
  
  #loading libraries
  require(foreach)
  require(doParallel)
  
  #create socket to save iteration
  # log.socket <- utils::make.socket(port=4000)
  # Log <- function(text, ...) {
  #   msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  #   cat(msg)
  #   write.socket(log.socket, msg)
  # }

  #NA string
  nas <- rep(NA, repetitions)
  
  #output dataframe
  out.df <- data.frame(
    method = nas,
    transformation = nas,
    distance = nas,
    diagonal = nas,
    paired.samples = nas,
    cases.to.modify = nas,
    remove.samples = nas,
    psi = nas,
    group = nas,
    stringsAsFactors = FALSE
  )
  
  #preparing cluster
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type="PSOCK")
  doParallel::registerDoParallel(my.cluster)
  
  #nrow params
  n.params <- nrow(params)
  
  #exporting cluster variables
  parallel::clusterExport(cl = my.cluster,
                          varlist = c('repetitionsPsi',
                                      'out.df',
                                      'params',
                                      'n.params',
                                      'repetitions'),
                          envir = environment())
  
  #iterating through repetitions
  out.df.list <- foreach::foreach(i = 1:n.params) %dopar% {
    
    # Log("Processing iteration %d of %d", i, n.params)
    
    #extracting params to environmnet
    list2env(x = as.list(params[i, ]), envir = environment())
    
    #filling out.df
    out.df[1:repetitions, "method"] <- method
    out.df[1:repetitions, "transformation"] <- transformation
    out.df[1:repetitions, "distance"] <- distance
    out.df[1:repetitions, "diagonal"] <- diagonal
    out.df[1:repetitions, "paired.samples"] <- paired.samples
    out.df[1:repetitions, "cases.to.modify"] <- cases.to.modify
    out.df[1:repetitions, "remove.samples"] <- remove.samples
    out.df[1:repetitions, "ignore.blocks"] <- ignore.blocks
    out.df[1:repetitions, "group"] <- i
    out.df[1:repetitions, "psi"] <- repetitionsPsi(
      x = data,
      transformation = transformation,
      distance = distance,
      diagonal = diagonal,
      paired.samples = paired.samples,
      cases.to.modify = cases.to.modify,
      remove.samples = remove.samples,
      ignore.blocks = ignore.blocks,
      repetitions = repetitions
    )
    
    return(out.df)
    
  }
  
  #stopping cluster
  parallel::stopCluster(my.cluster)
  
  #rbind
  out.df <- do.call("rbind", out.df.list)
  
  #aggregating by group and adding to params
  params$psi.mean <- aggregate(
    out.df[, c("psi", "group")], 
    by = list(out.df$group), 
    FUN = mean
  )$psi
  
  params$psi.sd <- aggregate(
    out.df[, c("psi", "group")], 
    by = list(out.df$group), 
    FUN = sd
  )$psi
  
  
  #generating a better name for method
  params[params$method == "elastic" & params$diagonal == TRUE & params$ignore.blocks == FALSE, "method"] <- "elastic - diagonal"
  params[params$method == "elastic" & params$diagonal == FALSE  & params$ignore.blocks == FALSE, "method"] <- "elastic - orthogonal"
  params[params$method == "elastic" & params$diagonal == TRUE & params$ignore.blocks == TRUE, "method"] <- "elastic - diagonal - no blocks"
  params[params$method == "elastic" & params$diagonal == FALSE  & params$ignore.blocks == TRUE, "method"] <- "elastic - orthogonal - no blocks"
  
  return(params)
  
}