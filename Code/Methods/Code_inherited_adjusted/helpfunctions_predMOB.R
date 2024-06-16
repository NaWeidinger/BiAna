### functions to assess permutation importance for a single tree
# function for lmtree
VI.lmtree <- function (object, newdat) 
{
  tids <- nodeids(object, terminal = TRUE)
  
  newdat$pnode <- predict(object, newdata = newdat, type = "node")
  
  objfuns <- unlist(sapply(1:nrow(newdat), FUN=function(i){
    xi <- newdat[i,]
    nodeid <- xi$pnode
    objfuns.i <- nodeapply(object, ids = nodeid, FUN = function(x) get_objfun_new( x, newdat=xi, var="pnode"))
    objfuns.i
  }))
  objfun.permvar <- lapply(attr(object$info$terms$partitioning, "term.labels"), FUN=function(var){
    newdat[, paste("pnode.perm", var, sep=".")] <- predictnode.perm(object, newdata = newdat, perm=var)
    objfuns <- unlist(sapply(1:nrow(newdat), FUN=function(i){
      xi <- newdat[i,]
      nodeid <- xi[paste("pnode.perm", var, sep=".")]
      objfuns.i <- nodeapply(object, ids = nodeid, FUN = function(x) get_objfun_new( x, newdat=xi, var=paste("pnode.perm", var, sep=".")))
      objfuns.i
    }))
    objfuns
  })
  y <- newdat[,getResponse(object$info$formula)]
  
  ret <- unlist(lapply(objfun.permvar, FUN=function(x) -mean(-(y-x)^2)+mean(-(y-objfuns)^2)))
  names(ret) <- attr(object$info$terms$partitioning, "term.labels")
  return(ret)
}

# function for glmtree
VI.glmtree <- function (object, newdat) 
{
  newdat <- data.frame(newdat)
  tids <- nodeids(object, terminal = TRUE)
  
  newdat$pnode <- predict(object, newdata = newdat, type = "node")
  
  objfuns <- unlist(nodeapply(object, ids = tids[sapply(tids, FUN=function(x) any(newdat$pnode==x))], 
                              FUN = function(x) get_objfun_new( x, newdat=newdat, var="pnode")))
  objfun.permvar <- lapply(attr(object$info$terms$partitioning, "term.labels"), FUN=function(var){
    newdat[, paste("pnode.perm", var, sep=".")] <- predictnode.perm(object, newdata = newdat, perm=var)
    objfuns <- unlist(nodeapply(object, ids = tids[sapply(tids, FUN=function(x) any(newdat[, paste("pnode.perm", var, sep=".")]==x))], 
                                FUN = function(x) get_objfun_new( x, newdat=newdat, var=paste("pnode.perm", var, sep="."))))
  })
  
  ret <- unlist(lapply(objfun.permvar, FUN=function(x) sum(objfuns)-sum(x)))
  names(ret) <- attr(object$info$terms$partitioning, "term.labels")
  return(ret)
}


# helpfunctions
getResponse <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response] 
}

predictnode.perm <- function(object, newdata = NULL, type = "node", perm=NULL) {
  predict.party(object, newdata = newdata, perm=perm)
}

get_objfun_new <- function(node, newdat, var) {
  nd <- subset(newdat, newdat[,var] == node$id)
  objfun(node$info$object, newdata=nd)
}

### functions to assess mean minimal depth for a list of trees
mindepth.mob.forest <- function (forest, mean_sample=c("all_trees", "top_trees", "relevant_trees")) {
  if(class(forest[[1]])[1]=="palmtree"){
    splitvarnames <- attr(forest[[1]]$tree$info$terms$partitioning,"term.labels")
  } else if (class(forest[[1]])[1] %in% c("glmtree", "lmtree", "modelparty")){
    splitvarnames <- attr(forest[[1]]$info$terms$partitioning,"term.labels")
  }
  
  importance_frame <- data.frame(variable = splitvarnames,
                                 stringsAsFactors = FALSE)
  ff <- function(x, results){
    if (depth(x)<1) {
      results <- rbind(results, data.frame(nodeID=x$id,
                                           left=0,
                                           right=0,
                                           splitvar=NA,
                                           splitpoint=NA))
    } else {
      results <- rbind(results, data.frame(nodeID=x$id,
                                           left=x$kids[[1]]$id,
                                           right=x$kids[[2]]$id,
                                           splitvar=character_split(split_node(x), data=forest[[1]]$data)$name,
                                           splitpoint=ifelse(is.null(x$split$breaks), NA, x$split$breaks)))
      if(depth(x$kids[[1]])>0){
        results <- ff(x$kids[[1]], results)
      } else {
        results <- rbind(results, data.frame(nodeID=x$kids[[1]]$id,
                                             left=0,
                                             right=0,
                                             splitvar=NA,
                                             splitpoint=NA))
      }
      if(depth(x$kids[[2]])>0){
        results <- ff(x$kids[[2]], results)
      } else {
        results <- rbind(results, data.frame(nodeID=x$kids[[2]]$id,
                                             left=0,
                                             right=0,
                                             splitvar=NA,
                                             splitpoint=NA))
      }
    } 
    return(results)
  }
  
  forest_table <- lapply(1:length(forest), FUN=function(i) {
    tree <- forest[[i]]
    results <- data.frame(nodeID=integer(),
                          left=integer(),
                          right=integer(),
                          splitvar=integer(),
                          splitpoint=double(), 
                          stringsAsFactors=FALSE)
    tree3 <- ff(node_party(tree), results)
    tree3 <- tree3[!duplicated(tree3$nodeID),]
    tree3 <- tree3[order(tree3$nodeID),]
    rownames(tree3) <- tree3$nodeID
    colnames(tree3) <- c("nodeID", "left daughter", "right daughter", "split var", "split point")
    if (!all(c("right daughter", "left daughter") %in% names(tree3))) {
      stop("The data frame has to contain columns called 'right daughter' and 'left daughter'!")
    }
    tree3$depth <- NA
    tree3$depth[1] <- 0
    if(nrow(tree3)>1){
      for (j in 2:nrow(tree3)) {
        tree3[j, "depth"] <- tree3[tree3[, "left daughter"] == 
                                     as.numeric(rownames(tree3[j, ])) | tree3[, "right daughter"] == 
                                     as.numeric(rownames(tree3[j, ])), "depth"] + 1
      }
    }
    nodesmat <- cbind(tree3, "tree"=i)
    return(nodesmat)
  })
  
  forest_table <- rbindlist(forest_table)
  forest_table <- forest_table[!is.na(forest_table$`split var`),]
  maxdepth <- dplyr::group_by(forest_table, tree) %>% dplyr::summarize(max(depth))
  min_depth_frame <- dplyr::group_by(forest_table, tree, `split var`)
  min_depth_frame <- dplyr::summarize(min_depth_frame, min(depth))
  colnames(min_depth_frame) <- c("tree", "variable", "minimal_depth")
  if(nrow(min_depth_frame)>0){
    importance_frame <- merge(importance_frame, randomForestExplainer:::measure_min_depth(min_depth_frame, 
                                                                                          mean_sample), by="variable", all = TRUE)
    importance_frame$mean_min_depth[is.na(importance_frame$mean_min_depth)] <- randomForestExplainer:::min_depth_count(min_depth_frame)[[3]]
  } else {
    importance_frame$mean_min_depth <- 1
  }
  return(importance_frame)
}  


