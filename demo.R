rm(list = ls())
source("functions.R")
source("getName.R")
library(bnlearn)
library(xlsx)
library(readxl)

pop_Size<-80;        #/* The number of Solutions*/ 种群个体个数
MaxAge <- 100;        # 最大感染年龄（50、100、300、500对应Sen5 -Sen8）
C0 <- 3;              # number of solutions have corona virus
max_Iter<-10;     #/*The number of cycles for foraging {a stopping criteria}*/最大迭代次数
orders <- 1;

file_name <- '100-1.xlsx'
set.seed(20240427)
gene_Num = 100
node_Name <- c(get_name("A", gene_Num), get_name("B", gene_Num), get_name("C", gene_Num))
xdata <- read_excel(file_name, sheet = 'Sheet1')
xdata <- xdata[, -1]
rownames(xdata) <- c(get_name("A", gene_Num), get_name("B", gene_Num), get_name("C", gene_Num))
xdata <- t(xdata)
xdata <- as.data.frame.array(xdata)
pear = read.xlsx(file_name, sheetIndex = 3)
rownames(pear) <- pear[,1]
pear = pear[,-1]
Lass = read.xlsx(file_name, sheetIndex = 4)
rownames(Lass) <- Lass[, 1]
Lass = Lass[, -1]
Gold <- read_excel(file_name, sheet = 'Gold')
Gold <- Gold[, -1]
rownames(Gold) <- get_name("G", gene_Num)

SpreadingRate <- 0.01;   # Spreading rate parameter传播速度参数（Sen1-Sen4对应0:005、0:05、0:01、0:5)，研究了基本繁殖率(BRr)对CHIO收敛性的影响）
runs <- 1;               #/*Algorithm can be run many times in order to see its robustness*/

swarm_score <- matrix(data = 0, nrow = 1, ncol = pop_Size);
Age <- matrix(data = -1, nrow = 1, ncol = pop_Size);
BestResults <- matrix(data = 0, nrow = runs, ncol = 1); # saving the best solution at each run

score_type = 'bic-g'
func_range <- c(1 :1)
for(i in func_range){
  if(i == 1){
    score_type = 'bic-g'
  }
  cat(sprintf("This is function %s", score_type), "\n")
  for(j in c(1: runs)){
    cat(sprintf("This is run %i", j), "\n")
    ###### initialization ######
    
    swarm = random_swarm(node_Name, pop_Size, gene_Num, 0.04)                 # (pop_Size) empty graphs
    swarm_score <- get_swarm_score(pop_Size, swarm, xdata, 'bic-g', pear, Lass, gene_Num, 2, 0.25)  # the score of initial graphs
    swarm_fitness <- get_swarm_fitness(pop_Size, swarm_score)           # gain fitness via score
    swarm_status <- matrix(data = 0, nrow = 1, ncol = pop_Size)         # initial status
    for (i in (1 : C0)){
      swarm_status[sample(1: pop_Size, 1, replace = TRUE)] = 1
    }
    ###### loop will begin immediately ######
    begin_time <- Sys.time()
    for(i in (1 : max_Iter)){
      # print(i)
      for(index in (1 : pop_Size)){
        NewSol <- swarm[[index]]
        is_Cornoa = 0
        confirmed = which(swarm_status == 1)
        normal = which(swarm_status == 0)
        immune = which(swarm_status == 2)
        for(order in (1 : orders)){
          rand <- runif(1)
          if ((rand < SpreadingRate / 3) && length(confirmed) > 0){
            ref <- sample(confirmed, 1, replace = TRUE)
            NewSol <- get_New_DBN(swarm[[index]], swarm[[ref]], node_Name, gene_Num)
            is_Cornoa = is_Cornoa + 1
          }else if((rand < SpreadingRate * 2 / 3) && length(normal) > 0){
            ref <- sample(normal, 1, replace = TRUE)
            NewSol <- get_New_DBN(swarm[[index]], swarm[[ref]], node_Name, gene_Num)
          }else if((rand < SpreadingRate) && length(immune) > 0){
            ref <- sample(immune, 1, replace = TRUE)
            NewSol <- get_New_DBN(swarm[[index]], swarm[[ref]], node_Name, gene_Num)
          }
        }
        #evaluate new solution
        scoreSol <- BIC_LP(NewSol, xdata, score_type, pear, Lass, gene_Num, 2, 0.25)
        fitnessSol <- get_swarm_fitness(1, scoreSol)
        # Update the curent solution  & Age of the current solution
        if (swarm_score[index]>scoreSol){
          swarm[[index]]=NewSol;
          swarm_fitness[index]=fitnessSol;
          swarm_score[index]=scoreSol;
        }
        else if(swarm_status[index]==1){
          Age[index] = Age[index] + 1;
        }
        # change the solution from normal to confirmed
        fitness_Mean <- sum(swarm_fitness) / length(swarm_fitness)
        if ((swarm_fitness[index] < fitness_Mean)&& swarm_status[index]==0 && is_Cornoa > 0){
          swarm_status[index] = 1;
          Age[index]=1;
        }
        # change the solution from confirmed to recovered
        if ((swarm_fitness[index] >= fitness_Mean)&& swarm_status[index]==1){
          swarm_status[index] <- 2;
          Age[index] <- 0;
        }
        # killed the current soluion and regenerated from scratch
        if(Age[index]>=MaxAge){
          swarm[[index]] <- empty.graph(node_Name, 1)
          swarm_status[index] <- 0;
          Age[index] <- -1;
        }
      }
    } 
    end_time <-Sys.time()
    print(end_time - begin_time)
    ind = which(swarm_score == min(swarm_score))
    score_result <- BIC_LP(swarm[[ind]], xdata, score_type, pear, Lass, gene_Num, 2, 0.25)
    print("Score_result")
    print(score_result)
    resmat <- get_res_mat(swarm[[ind]], gene_Num)
    rres = calResult(resmat, Gold, gene_Num)
    print(rres)
  }
}
