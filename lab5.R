library(igraph)
library(igraphdata)
library(clustAnalytics)
data(karate,package="igraphdata")
wc <- walktrap.community(karate)
modularity(wc)
#[1] 0.3532216
unname(membership(wc)) ## try membership(wc)
#[1] 1 1 2 1 5 5 5 1 2 2 5 1 1 2 3 3 5 1 3 1 3 1 3 4 4 4 3 4 2 3 2 2 3 3
plot(wc, karate)

data(karate,package="igraphdata")
fc <- fastgreedy.community(karate)
dendPlot(fc)

as_adjacency_matrix(as.undirected(karate,mode = "each"))

evaluate_significance(karate,alg_list=list(Louvain=cluster_louvain,
                                           "label prop"= cluster_label_prop,
                                           walktrap=cluster_walktrap),
                      gt_clustering=V(karate)$Faction)

B <- matrix(c(1, 0.2, 0.2, 1), ncol=2)
G <- barabasi_albert_blocks(m=4, p=c(0.5, 0.5), B=B, t_max=100,
                            type="Hajek", sample_with_replacement = FALSE)
plot(G, vertex.color=(V(G)$label),vertex.label=NA,vertex.size=10)


#intersection: count all elements of A that are in B

#union: sum of cardinalities - intersection

jaccard_idx <- function(A, B){
  
  int <- c()
  for(elem in A){
    if(elem %in% B){
      int <- append(elem,int) 
      }
  }
  int <- length(int)
  union = length(A) + length(B) - int
  jaccard_idx = int/union
  return(jaccard_idx)
  
}

#jaccard_sim

jaccard_sim <- function(Ca, Cb){ #Ca on rows, Cb on columns
  tab <- matrix(c(1:(length(Ca)*length(Cb))), ncol = length(Cb), byrow=TRUE)
  colnames(tab) = c(1:length(Cb))
  rownames(tab) <- c(1:length(Ca))
  for (i in 1:(length(Cb))){
    for(j in 1:(length(Ca))){ #Since matrix go by columns first, loops are inverted
      tab[j+length(Ca)*(i-1)] = jaccard_idx(unlist(Ca[j]),unlist(Cb[i])) #j+length(Ca)*(i-1) is the index calculation
      

    }
    #print(tab)
  }
  return(tab)
}

match_clusters <- function(JS, name1, name2){
  
  cols <- ncol(JS) 
  rows <- nrow(JS)
  
  match <- c()
  lab <- c()
  
  for(i in 1:cols){
    s <- 0
    idxi <- 0
    idxj <- 0
    for(j in 1:rows){
      if(s < JS[j + rows*(i-1)]){
        s = JS[j + rows*(i-1)]
        idxi <- i
        idxj <- j
      }
    }
    match <- append(match,s)
    lab <- append(lab,paste("(",name1,".",idxj,name2,".",idxj,")"))
  }
  m <- matrix(match, ncol = length(match), byrow=TRUE)
  colnames(m) = lab
  return(m)
}


fc <- fastgreedy.community(karate)
wc <- walktrap.community(karate)
JS <- jaccard_sim(fc,wc)
print(JS)
ncol(JS) 
nrow(JS)
match_clusters(JS,"FC","WC")
