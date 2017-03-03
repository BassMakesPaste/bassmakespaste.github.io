## Build essential functions
read_table_shiny = function(inFile, header=FALSE)
{
  out = NULL
  if (!is.null(inFile$datapath))
  {
    if(file_ext(inFile$name) != 'csv')
    {
      if(is.null(inFile))
        return(NULL)
      
      ext = file_ext(inFile$name)
      
      if (file.exists(inFile$datapath))
        file.rename(inFile$datapath, paste(inFile$datapath, ext, sep="."))
      
      if (ext == 'xlsx')
        out = read_excel(paste(inFile$datapath, ext, sep="."), col_names=header)
      else
        out = read.csv(inFile$datapath, header=header)
    } else {
      out = read.csv(inFile$datapath, header=header)
    }
  }
  return(out)
}

permn <- function(x) {
    if (length(x) == 1) {
        return(x)
    } else {
        r = matrix(nrow = 0, ncol = length(x))
        for (i in seq_along(x)) {
            r = rbind(r, cbind(x[i], Recall(x[-i])))
        }
        return(r)
    }
}


## Imports
library(shiny)
library(readxl)
library(tools)

shinyServer(function(input, output) {
  output$out <- renderPlot({
    if (!is.null(input$complementation$datapath))
    {
        complementation = as.matrix(read_table_shiny(input$complementation, header=FALSE))
        checkMap = input$checkMap
        if (!is.null(input$mapping$datapath))
            mapping = as.matrix(read_table_shiny(input$mapping, header=FALSE))

        ######################################################################################################
        # BEGIN CORE FUNCTION
        {
            # Preprocess matrices to make sure that they're goodo.
            {
                # Check to make sure complementation and mapping tables are square.
                if (nrow(complementation) != ncol(complementation))
                {
                    stop('Matrices are not square: aborting.')
                }
                if (nrow(mapping) != ncol(mapping))
                {
                    stop('Matrices are not square: aborting.')
                }
                if (nrow(mapping) != ncol(complementation))
                {
                    stop('Matrices are not of the same dimensions: aborting.')
                }
                
                nmutants = nrow(complementation)
                # Check if they are lower/upper triange, then convert them if so.
                # Store 4 binary states as one variable by multiplying them by factors of 2.
                status = ((mapping[1,nmutants] == '') + 2*(mapping[nmutants,1] == '') + 4*(mapping[nmutants,nmutants] == '') + 8*(mapping[1,1] == ''))
                #Matrix type:  Upper triangle                 Upper triangle                  Lower right triangle               Upper right triangle
               
                # Keep in mind that right triangles are thorougly untested--flip your triangles if you're going to use them!
                if (status %in% c(0,1,2,4,8))
                {
                    # Depending on the matrix type, we pick a function to operate on it.
                    f = switch(paste('p',as.character(status),sep=''),
                                p0 = function(x) return(x),
                                p1 = function(x) {v = x[lower.tri(x)]; x = t(x); x[lower.tri(x)] = v; return(x)},
                                p2 = function(x) {v = x[upper.tri(x)]; x = t(x); x[upper.tri(x)] = v; return(x)},
                                p4 = function(x) {x = x[c(nmutants:1),,drop = FALSE]; v = x[lower.tri(x)]; x = t(x); x[lower.tri(x)] = v; x = x[c(nmutants:1),,drop = FALSE]; return(x)},
                                p8 = function(x) {x = x[,c(nmutants:1),drop = FALSE]; v = x[upper.tri(x)]; x = t(x); x[upper.tri(x)] = v; x = x[c(nmutants:1),,drop = FALSE]; return(x)}
                          )
                    mapping = f(mapping)
                } else {
                    stop('Multiple matrix corners are empty. Please use a lower/upper triangle or square matrix.')
                }
                
                status = ((complementation[nmutants,1] == '') + 2*(complementation[1,nmutants] == '') + 4*(complementation[nmutants,nmutants] == '') + 8*(complementation[1,1] == ''))
                #Matrix type:  Lower triangle                 Upper triangle                  Upper right triangle
                    # Keep in mind that right triangles are thorougly untested--flip your triangles if you're going to use them!
                if (status %in% c(0,1,2,4,8))
                {
                    # Depending on the matrix type, we pick a function to operate on it.
                    f = switch(paste('p',as.character(status),sep=''),
                                p0 = function(x) return(x),
                                p1 = function(x) {v = x[upper.tri(x)]; x = t(x); x[upper.tri(x)] = v; return(x)},
                                p2 = function(x) {v = x[lower.tri(x)]; x = t(x); x[lower.tri(x)] = v; return(x)},
                                p4 = function(x) {x = x[,c(nmutants:1),drop = FALSE]; v = x[upper.tri(x)]; x = t(x); x[upper.tri(x)] = v; x = x[c(nmutants:1),,drop = FALSE]; return(x)},
                                p8 = function(x) {x = x[c(nmutants:1),,drop = FALSE]; v = x[lower.tri(x)]; x = t(x); x[lower.tri(x)] = v; x = x[c(nmutants:1),,drop = FALSE]; return(x)}
                          )
                    complementation = f(complementation)
                } else {
                    stop('Multiple matrix corners are empty. Please use a lower/upper triangle or square matrix.')
                }
                
                # Change +/- to 0/1, if needed.
                complementation[complementation == '-'] = '0'
                complementation[complementation == '+'] = '1'
                complementation = matrix(as.numeric(complementation),nrow=nmutants,ncol=nmutants) 
                
                mapping[mapping == '-'] = '0'
                mapping[mapping == '+'] = '1'
                mapping = matrix(as.numeric(mapping),nrow=nmutants,ncol=nmutants)
                
                ##if (!isSymmetric(mapping))
                ##{
                ##    stop('Matrices are not symmetric: aborting.')
                ##}
                if (!isSymmetric(complementation))
                {
                    stop('Matrices are not symmetric: aborting.')
                }
            }    

            # Establish some vars.
            geneRegions = list()
            pointMutants = vector()
            nonpointMutants = vector()
            nonpointOccurrence = vector()

            # Analyze complementation table and find gene regions, as well as point mutants.
            {    
                ## Get information from mapping and complementation tables.
                for(i in 1:nmutants)
                {
                    assign(paste('anti', i, sep=''), which(complementation[i,] == 0))
                }
                
                for (i in 1:nmutants)
                {
                    antiCurrent = get(paste('anti',i,sep=''))
                    isGeneRegion = TRUE
                    for (j in antiCurrent)
                    {
                        if (length(get(paste('anti',j,sep=''))) < length(antiCurrent))
                        {
                            isGeneRegion = FALSE
                        }
                    }
                    if (isGeneRegion == TRUE)
                    {
                        geneRegions[[i]] = antiCurrent
                        pointMutants = c(pointMutants,i)
                    } else {
                        geneRegions[[i]] = 0
                        nonpointMutants = c(nonpointMutants,i)
                    }
                }
                geneRegions = unique(geneRegions)
            }

            # Go through the mapping table to weed out single-gene deletion mutants.
            if (checkMap) {
                for(i in 1:length(pointMutants))
                {
                    assign(paste('maps', pointMutants[i], sep=''), which(mapping[pointMutants[i],] == 1))
                    mapCurrent = which(mapping[pointMutants[i],] == 1)
                    mapCurrent = mapCurrent[mapCurrent %in% pointMutants]
                    antiCurrent = get(paste('anti',pointMutants[i],sep=''))
                    antiCurrent = antiCurrent[antiCurrent %in% pointMutants]
                    antiCurrent = antiCurrent[antiCurrent != pointMutants[i]]
                    
                    pointMutation = FALSE
                    if (length(antiCurrent) > 0)
                        {
                        for(j in 1:length(antiCurrent))
                        {
                            if(antiCurrent[j] %in% mapCurrent)
                            {
                                pointMutation = TRUE
                            }
                        }
                        
                        if (pointMutation == FALSE)
                        {
                            nonpointMutants = c(nonpointMutants, pointMutants[i])
                            pointMutants = pointMutants[pointMutants != pointMutants[i]]
                            nonpointOccurrence = c(nonpointOccurrence,1)
                        }
                    }
                }
            }

            # Generate object which we use to order point mutants.
            overlapRegions = list()
            for (i in 1:length(geneRegions))
            {
                if (!identical(geneRegions[[i]],0))
                {
                   overlapRegions[[length(overlapRegions)+1]] = list(pointMutants=geneRegions[[i]][geneRegions[[i]] %in% pointMutants],nonpointMutants=geneRegions[[i]][geneRegions[[i]] %in% nonpointMutants])
                }
            }

            for (i in 1:length(nonpointMutants))
            {
                nonpointOccurrence[i] = sum(unlist(overlapRegions) == nonpointMutants[i])
            }

            # Get all possible perms, then iterate until we can order the chromosome in a 
            # way which satisfies our tables.
            # I'm too lazy to figure out how to do this properly.
            unsolved = TRUE
            nums = permn(c(1:length(overlapRegions)))
            cc = 1
            while(unsolved == TRUE && cc <= nrow(nums))
            {
                sampleRegions = list()
                for (i in 1:ncol(nums))
                {
                    sampleRegions[[i]] = overlapRegions[[nums[cc,i]]]
                }
                
                hasbeenused = vector()
                used = vector()
                endusage = vector()
                solver = vector()
                
                for (i in 1:length(overlapRegions))
                {
                    used = sampleRegions[[i]]$nonpointMutants
                    
                    for (j in 1:length(used))
                    {
                        if ((used[j] %in% hasbeenused) && (used[j] %in% endusage))
                        {
                            solver = c(solver,TRUE)
                        } else {
                            solver = c(solver,FALSE)
                        }
                    }
                    hasbeenused = unique(c(used,hasbeenused))
                    endusage = c(endusage,hasbeenused[!(unlist(lapply(hasbeenused, function(x) return(x %in% used))))])
                }
                unsolved = as.logical(sum(solver))
                cc = cc + 1
            }

            # Plot our sampleRegions.

            par(mar=c(0,0,0,0))
            plot(0, type = 'n' , xaxt = 'n', yaxt = 'n', xlab = FALSE, ylab = FALSE, ann = FALSE, xlim=c(0,120), ylim=c(0,100))
            lines(c(0,120),c(50,50))
            boxes = lapply(c(1:nmutants),function(x) return(NULL))
            for (i in 1:length(sampleRegions))
            {
                y = 50
                x = ((100*i)/length(sampleRegions))
                p = sampleRegions[[i]]$pointMutants
                text(x=x,y=y, paste(replicate(length(p),'X'),collapse=' '),offset=0)
                text(x=x,y=y+10, paste(p,collapse=' '),offset=0)
                
                n = sampleRegions[[i]]$nonpointMutants
                for (j in n)
                {
                    if (is.null(boxes[[j]]))
                    {
                        boxes[[j]]$start = i
                    }
                    for(k in i:length(sampleRegions))
                    {
                        if (!(j %in% sampleRegions[[k]]$nonpointMutants) && is.null(boxes[[j]]$end))
                        {
                            boxes[[j]]$end = k-1
                        }
                       
                    }
                    if (is.null(boxes[[j]]$end))
                    {
                        boxes[[j]]$end = k
                    }
                }
            }

            for (i in 1:length(nonpointMutants))
            {
                x1 = (100*boxes[[nonpointMutants[i]]]$start)/length(sampleRegions)
                y1 = 44 - ((i-1) * 10)
                x2 = (100*boxes[[nonpointMutants[i]]]$end)/length(sampleRegions)
                y2 = y1 - 8
                
                rect(x1,y2,x2,y1)
                text(x=x1+((x2-x1)/2),y=y2+((y1-y2)/2), nonpointMutants[i])
            }
            
            
        }    
        ######################################################################################################
        # END CORE FUNCTION
        
        
    } else {
      stop('Please upload a table. Accepted formats are: csv, txt, xlsx\n')
    }
    
  })
})
