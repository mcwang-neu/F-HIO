checkCircleInMat <- function(graph, nNode){
  
  nParentNum = getParentsNum(graph, nNode)
  NodeStack = matrix(data = -1, nrow = 1, ncol = nNode)
  NodeVisit = matrix(data = 0, nrow = 1, ncol = nNode)
  top = 0
  i = 1
  while(i <= nNode){
    if(nParentNum[i] == 0){
      top = top + 1
      NodeVisit[i] = 1
      NodeStack[top] = i
    }
    i = i + 1
  }
  while(top != 0){
    flag = NodeStack[top]
    top = top - 1
    NodeVisit[flag] = 1
    
    j = 1
    while(j <= nNode){
      if(graph[flag, j] == 1){
        nParentNum[j] = nParentNum[j] - 1
        if(nParentNum[j] == 0 && NodeVisit[j] == 0){
          top = top + 1 
          NodeStack[top] = j
        }
      }
      j = j + 1
    }
  }
  i = 1
  while(i <= nNode){
    if(NodeVisit[i] == 0){
      return(TRUE)
    }
    i = i + 1
  }
  return(FALSE)
}

getParentsNum <- function(graphs, nNode){
  ret = matrix(data = 0, nrow = 1, ncol = nNode)
  j = 1
  while(j <= nNode){
    i = 1
    while(i <= nNode){
      if(graphs[i,j] == 1){
        ret[j] = ret[j] + 1
      }
      i = i + 1
    }
    j = j + 1
  }
  return(ret)
}


BIC_LP<-function(net, xdata, score_type, pear, lasso, vertexNum, m, gamma){
  sum = score(net, xdata, type = score_type)
  if(m == 0){
    return(sum)
  }
  secondMat = net2SecondMat(net, vertexNum)
  first_Mat = net2First_Mat(net, vertexNum)
  for (i in 1:vertexNum) {
    for(j in 1 : vertexNum){
      if(i == j){
        next
      }
      if(secondMat[i,j] == 1){
        sum = sum + log(m * lasso[i,j] + gamma)
      }else{
        sum = sum + log(m * (1 - lasso[i,j]) + gamma)
      }
      if(first_Mat[i,j] == 1){
        sum = sum + log(m * lasso[i + vertexNum, j] + gamma)
      }else{
        sum = sum + log(m * (1 - lasso[i + vertexNum, j]) + gamma)
      }
    }
  }
  # print(sum)
  return(-sum)
}


get_swarm_score <- function(pop_Size, swarm, xdata, score_type, pear, lasso, vertexNum, m, gamma){
  res = c()
  for(i2 in (1: pop_Size)){
    this_score <- BIC_LP(swarm[[i2]], xdata, score_type, pear, lasso, vertexNum, m, gamma)
    res = c(res, this_score)
    # print(this_score)
  }
  #print(res)
  return(res)
}


get_swarm_fitness <- function(pop_Size, swarm_score){
  res = c()
  this_fitness <- NULL
  for (i in (1 :  pop_Size)){
    # print(i)
    if(swarm_score[i] > 0){
      this_fitness <- 1 / (swarm_score[i] + 1)
    }else{
      this_fitness <- 1 + abs(swarm_score[i])
    }
    res = c(res, this_fitness)
  }
  return(res)
}


net2StaticMat <- function(net, vertexNum){
  res = matrix(0, nrow = vertexNum, ncol = vertexNum)
  len = length(net[["arcs"]]) / 2
  
  i = 1
  while(i <= len){
    from = net[["arcs"]][i]
    to = net[["arcs"]][i + len]
    
    nfrom = nchar(from)
    nto = nchar(to)
    
    if(substr(from, 1, 1) == "C"){
      from = substr(from, 2, nfrom)
      to = substr(to, 2, nto)
      index1 = which(c(1:vertexNum) == from)
      index2 = which(c(1:vertexNum) == to)
      res[index1, index2] = 1
    }
    i = i + 1
  }
  
  return(res)
}


net2First_Mat <- function(net, vertexNum){
  res = matrix(0, nrow = vertexNum, ncol = vertexNum)
  len = length(net[["arcs"]]) / 2
  i = 1
  while (i <= len){
    from = net[["arcs"]][i]
    to = net[["arcs"]][i + len]
    nfrom = nchar(from)
    nto = nchar(to)
    if(substr(from, 1, 1) == "B"){
      from = substr(from, 2, nfrom)
      to = substr(to, 2, nto)
      index1 = which(c(1:vertexNum) == from)
      index2 = which(c(1:vertexNum) == to)
      res[index1, index2] = 1
    }
    i = i + 1
  }
  return(res)
}

net2SecondMat <- function(net, vertexNum){
  res = matrix(0, nrow = vertexNum, ncol = vertexNum)
  len = length(net[["arcs"]]) / 2
  i = 1
  while(i <= len){
    from = net[["arcs"]][i]
    to = net[["arcs"]][i + len]
    
    nfrom = nchar(from)
    nto = nchar(to)
    
    if(substr(from, 1, 1) == "A"){
      from = substr(from, 2, nfrom)
      to = substr(to, 2, nto)
      index1 = which(c(1:vertexNum) == from)
      index2 = which(c(1:vertexNum) == to)
      res[index1, index2] = 1
    }
    i = i + 1
  }
  return(res)
}

Mat2DBN <- function(first_Mat, second_Mat, staticMat, nodeName, vertexNum){
  ret = empty.graph(nodeName)
  for(i in (1 : vertexNum)){
    for(j in (1 : vertexNum)){
      if(first_Mat[i, j] == 1){
        ret = set.arc(ret, nodeName[vertexNum + i], nodeName[vertexNum * 2 + j])
      }
      if(second_Mat[i,j] == 1){
        ret = set.arc(ret, nodeName[i], nodeName[vertexNum * 2 + j])
      }
      if(staticMat[i,j] == 1){
        ret = set.arc(ret, nodeName[vertexNum * 2 + i], nodeName[vertexNum * 2 + j])
      }
    }
  }
  return(ret)
}

get_New_DBN <- function(previous, reference, node_Name, vertexNum){
  previous_mat_first <- net2First_Mat(previous, vertexNum)
  previous_mat_second<- net2SecondMat(previous, vertexNum)
  previous_mat_static<- net2StaticMat(previous, vertexNum)
  previous_mat_static_2 <- net2StaticMat(previous, vertexNum)
  reference_mat_first <- net2First_Mat(reference, vertexNum)
  reference_mat_second<- net2SecondMat(reference, vertexNum)
  reference_mat_static<- net2StaticMat(reference, vertexNum)
  alpha <- 0.9                    #传染因子
  beta <- 0.05                    #突变因子
  for (i in (1 : vertexNum)){
    for(j in (1 : vertexNum)){
      if(i == j){
        next
      }
      rand_1 <- runif(1)
      rand_2 <- runif(1)
      rand_0 <- runif(1)
      if(rand_1 > alpha){
        previous_mat_first[i,j] <- reference_mat_first[i,j]
      }
      else if(rand_2 < beta){
        # print("发生突变")
        previous_mat_first[i,j] <- 1 - previous_mat_first[i,j]
      }
      if(rand_2 > alpha){
        previous_mat_second[i,j] <- reference_mat_second[i,j]
      }
      else if(rand_2 < beta){
        # print("发生突变")
        previous_mat_second[i,j] <- 1 - previous_mat_second[i,j]
      }
      if(rand_0 > alpha){
        previous_mat_static[i,j] <- reference_mat_static[i,j]
      }
      else if(rand_0 < beta){
        # print("发生突变")
        previous_mat_static[i,j] <- 1 - previous_mat_static[i,j]
      }
    }
  }
  if (checkCircleInMat(previous_mat_static, vertexNum)){
    return(Mat2DBN(previous_mat_first, previous_mat_second, previous_mat_static_2, node_Name, vertexNum))
  }
  return(Mat2DBN(previous_mat_first, previous_mat_second, previous_mat_static, node_Name, vertexNum))
}


random_swarm <- function(node_Name, pop_Size, vertexNum, alpha){# alpha 代表系数程度 0.05
  swarm <- empty.graph(node_Name, pop_Size)
  for(index in (1 : pop_Size)){
    arc_num <- round(alpha * vertexNum * vertexNum)
    while(arc_num > 0){
      a <- sample(1 : vertexNum, 2, replace = FALSE)
      swarm[[index]] = set.arc(swarm[[index]], node_Name[a[1]], node_Name[vertexNum * 2 + a[2]])
      arc_num = arc_num - 1
    }
  }
  return(swarm)
}

get_res_mat <- function(resDBN, vertexNum){
  res <- matrix(data = 0, nrow = vertexNum, ncol = vertexNum)
  mat1 <- net2First_Mat(resDBN, vertexNum)
  mat2 <- net2SecondMat(resDBN, vertexNum)
  for(i in (1 : vertexNum)){
    for(j in (1 : vertexNum)){
      if(mat1[i,j] ==1 || mat2[i,j] == 1){
        res[i,j] =1
      }
    }
  }
  return(res)
}

calResult <- function(resMat, Goldmat, vertexNum){
  
  # 和金标准(Goldmat)进行比较
  # resMat为片内和片间网络组合的二维表
  # Goldmat为金标准网络
  # vertexNum 为基因个数
  
  TP = 0
  TN = 0
  FP = 0
  FN = 0
  for(i in 1 : vertexNum){
    for(j in 1 :vertexNum){
      if(i == j){
        next
      }
      if(resMat[i,j] == 1 && Goldmat[i,j] == 1){
        TP = TP + 1
      }else if(resMat[i,j] == 1 && Goldmat[i,j] == 0){
        FP = FP + 1
      }else if(resMat[i,j] == 0 && Goldmat[i,j] == 0){
        TN = TN + 1
      }else{
        FN = FN + 1
      }
    }
  }
  
  Pre0 = round(TP/(TP+FP), 4)
  Recall0 = round(TP/(TP+FN), 4)
  Acc0 = round((TP+TN)/(TP+TN+FN+FP), 4)
  FScore0 = round(2 * Pre0 * Recall0 / (Pre0 + Recall0), 4)
  BAC = round((TP/(TP + FN) + TN/(TN + FP)) / 2, 4)
  res = c(TP, TN, FP, FN, Pre0, Recall0, Acc0, FScore0, BAC)
  mat0 = matrix(data = res, nrow = 1, ncol = 9, byrow = FALSE,dimnames = NULL)
  
  return(mat0)
}


