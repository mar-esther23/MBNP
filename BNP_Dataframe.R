#
attractorToDataframe <- function(attr, sep="/", node.names=NULL, Boolean=FALSE) {
    if (Boolean) {
          if (is(attr, "AttractorInfo")) {
        if (is.null(node.names)) node.names <- attr$stateInfo$genes
        attr <- sapply(attr$attractors, function(a) paste(a$involvedStates, collapse=sep) )
      }
      if (is.null(node.names)) stop("Invalid node.names")
      if (is.list(attr)) { attr <- unlist(attr) }

      df <- data.frame(matrix(ncol=length(node.names)+2, nrow=0))
      for (i in seq(length(attr))) {
        s <- attr[i]
        if(is.character(s)) s<-unlist(strsplit(s,sep))
        for (j in seq(length(s))) {
          df <- rbind(df, c(attractor=i,state=j,int2binState(s[j],node.names)))
        }
      }
      colnames(df)<-c('attractor','state',node.names)
      return(df)
    }
    
    else {
                if (!is(attr, "AttractorInfo")) { stop("Error: non-valid attractor") }
        attr <- attr$attractors
                attr.properties <- vector("list", length(attr[[(1)]]))
        names(attr.properties) <- names(attr[[(1)]])
        
        for (n in names(attr.properties) ) {             attr.properties[[n]] <- lapply(attr, function(a) a[[n]]) 
                        ncol <- max(sapply(attr.properties[[n]], length))
            if ( ncol > 1) {                   attr.properties[[n]] <- sapply(attr.properties[[n]], function(a) {
                  paste(as.character(a), collapse=sep)
                  })}
            attr.properties[[n]] <- unlist(attr.properties[[n]])
                    }
        return(data.frame(attr.properties, stringsAsFactors=F))
    }
}



attractorListToDataframe <- function(attr.list, sep='/', returnDataFrame=c('occurrence','basinSize'), ...) {
      returnDataFrame <- match.arg(returnDataFrame)
    attr.list <- lapply(attr.list, attractorToDataframe, sep)
    if ( is.null(names(attr.list)) ) names(attr.list) <- 1:length(attr.list)
    for (n in names(attr.list)) {
    rownames(attr.list[[n]]) <- attr.list[[n]]$involvedStates     if (returnDataFrame=='occurrence') attr.list[[n]] <- replace(attr.list[[n]],!is.na(attr.list[[n]]),1)
    attr.list[[n]]$involvedStates <- NULL     names(attr.list[[n]]) <- n       } 
  
    attr.df <- Reduce(function(x, y){
    df <- merge(x, y, by= "row.names", all=TRUE)
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
    return(df)
  }, attr.list)
  attr.df
}

dataframeToAttractor <- function(df, fixedGenes) {
  bin2intState <- function(x){ 
    x <- rev(x)
    sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
  }
  attractors = vector("list", length = max(df$attractor))
    for(i in 1:nrow(df)) {
    row <- df[i,]
    n <- row[[1]]
    state <- bin2intState(row[c(-1,-2)])
    attractors[[ n ]] <- c(attractors[[ n ]], state)
  }
  for (i in seq(length(attractors))) {
    l = length(attractors[[i]])
    attractors[[i]] <- list(
      involvedStates = array(attractors[[i]], dim=c(1,l)),
      basinSize = NA
    )  
  }
  node.names <- colnames(df)[c(-1,-2)]
  if (missing(fixedGenes)) {
    fixedGenes <- rep(-1, length(attr$stateInfo$genes))
    names(fixedGenes) <- node.names
  }
  stateInfo = list( genes = node.names, fixedGenes = fixedGenes )
  
  result <- list( stateInfo = stateInfo,attractors = attractors )
  class(result) <- "AttractorInfo"
  result
}
      
      
      

aggregateByLabel <- function(df, node.names, label.rules, sep='/') {
  labels <- lapply(rownames(df), function(state) {
    state <- as.numeric(unlist(strsplit(state, sep)))
    label <- lapply(state, function(s) {
      l <- labelState(s, node.names, label.rules)
    })  
    label <- paste(label, collapse=sep)
  })
  labels <- unlist(labels)
  df[is.na(df)] <- 0
  df <- aggregate(df, list(labels), sum)
  rownames(df) <- df$Group.1
  df$Group.1 <- NULL
  as.data.frame(df)    
}

      
countChangesDataframe <- function(df, reference='WT', axis=2) {
    df <- df/df
    df[is.na(df)] <- 0
    df <- df-df[[reference]]
    new <- apply(df, axis, function(x) sum(x==1))
    lost <- apply(df, axis, function(x) sum(x==-1))
    rbind(new, lost)
}
