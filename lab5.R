library(igraph)
library(igraphdata)
library(clustAnalytics)ç
library(xtable)
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

match_clusters <- function(JS, name1, name2){ #Matching of the most similar cluster labeling according to the jaccard index
  
  cols <- ncol(JS) 
  rows <- nrow(JS)
  
  match <- c()
  lab <- c()
  
  for(i in 1:rows){
    s <- 0
    idxi <- 0
    idxj <- 0
    for(j in 1:cols){
      print(JS[i + rows*(j-1)])
      if(s < JS[i + rows*(j-1)]){
        s = JS[i + rows*(j-1)]
        idxi <- i
        idxj <- j
        print(j)
      }
    }
    match <- append(match,s)
    lab <- append(lab,paste("(",name1,".",idxi,name2,".",idxj,")"))
  }
  m <- matrix(match, ncol = length(match), byrow=TRUE)
  colnames(m) = lab
  return(m)
}

Wmean <- function(MC, GT){ #Mean of the match_clusters vector taking into account the ground truth cluster size
  s <- 0
  total <-0
  for(i in 1:length(GT)){
    total <- total + length(unlist(GT[i]))
  }
  
  for(i in 1:length(MC)){
    s <- s + MC[i]/(total/length(unlist(GT[i])))
  }
  
  return(s)
}

make_gt <- function(gt, V){
  l <- list()
  for(i in 1:max(gt)){
      lc <- list() 
      for(j in 1:length(gt)){
        #cat("gt_j: ", gt[j], "i: ",i)
        if(gt[j]==i){
          
          lc <- append(lc, V[j]$name)
          print(lc)
        }
      }  
     
      print(lc)
      l <- append(l, list(lc))
  }
  print(l)
  for(i in 1:length(l)){
    l[i] = list(unlist(l[i]))
  }
  #print(l)  
  return(l)
  
  
}
  
make_gt_numeric <- function(gt){
  l <- list()
  for(i in 1:max(gt)){
    lc <- list() 
    for(j in 1:length(gt)){
      #cat("gt_j: ", gt[j], "i: ",i)
      if(gt[j]==i){
        
        lc <- append(lc, j)
        print(lc)
      }
    }  
    
    print(lc)
    l <- append(l, list(lc))
  }
  print(l)
  for(i in 1:length(l)){
    l[i] = list(unlist(l[i]))
  }
  #print(l)  
  return(l)
  
  
} 
  


# fc <- fastgreedy.community(karate)
# wc <- walktrap.community(karate)
# wc
# fc
# unname(membership(wc))
# V(karate)$Faction
# JS <- jaccard_sim(fc,wc)
# print(JS)
# ncol(JS) 
# nrow(JS)
# MC <- match_clusters(JS,"FC","WC")
# Wmean(MC,fc)

#Evaluation

#Karate

data(karate,package="igraphdata")

gt <- V(karate)$Faction
gt
gtc <- make_gt(gt,V(karate))


louc <- cluster_louvain(karate)
louc

JS <- jaccard_sim(gtc,louc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lou")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)


labc <- cluster_label_prop(karate)

JS <- jaccard_sim(gtc,labc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lab")
MC
print(xtable(MC, type = "latex"))
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)


wc <- walktrap.community(karate)
wc

JS <- jaccard_sim(gtc,wc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","walk")
MC
print(xtable(MC, type = "latex"))

Wmean(MC,gtc)


ebc <- cluster_edge_betweenness(karate)
ebc

JS <- jaccard_sim(gtc,ebc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","EB")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)
       
imc <- cluster_infomap(karate)
imc

JS <- jaccard_sim(gtc,imc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","EB")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

#Barabasi albert

B <- matrix(c( 1, 0.15, 0.1, 0.2,
               0.15, 1, 0.1, 0.2,
               0.15, 0.1, 1, 0.2,
               0.2, 0.15, 0.1, 1), ncol=4)
G <- barabasi_albert_blocks(m=4, p=c(0.1, 0.3, 0.4, 0.2), B=B, t_max=200,
                            type="Hajek", sample_with_replacement = FALSE)
plot(G, vertex.color=(V(G)$label),vertex.label=NA,vertex.size=10)

tab <-evaluate_significance(G,alg_list=list(Louvain=cluster_louvain,
                                                         "label prop"= cluster_label_prop,
                                                         walktrap=cluster_walktrap, 
                                                         "Edge bet." = cluster_edge_betweenness,
                                                         infomap=cluster_infomap)
                                                      ,gt_clustering=V(G)$label)

print(xtable(tab, type = "latex"))


gt <- V(G)$label
gt
V <- V(G)
gtc <- make_gt_numeric(gt)
length(gtc)

louc <- cluster_louvain(G)
louc

JS <- jaccard_sim(gtc,louc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lou")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)


labc <- cluster_label_prop(G)
labc

JS <- jaccard_sim(gtc,labc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lab")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

wc <- walktrap.community(G)
wc

JS <- jaccard_sim(gtc,wc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","walk")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

ebc <- cluster_edge_betweenness(G)
ebc

JS <- jaccard_sim(gtc,ebc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","EB")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

imc <- cluster_infomap(G)
imc

JS <- jaccard_sim(gtc,imc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","EB")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)


#ENRON

data(enron,package="igraphdata")

# evaluate_significance(as.undirected(enron),alg_list=list(Louvain=cluster_louvain,
#                                             "label prop"= cluster_label_prop,
#                                             walktrap=cluster_walktrap, "Edge bet." = cluster_edge_betweenness,
#                                           infomap=cluster_infomap))
gtc <- cluster_infomap(G)
gtc

louc <- cluster_louvain(G)
louc
JS <- jaccard_sim(gtc,louc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lou")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)


labc <- cluster_label_prop(G)
labc
JS <- jaccard_sim(gtc,labc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lab")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)


wc <- walktrap.community(G)
wc
JS <- jaccard_sim(gtc,wc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Walk")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

ebc <- cluster_edge_betweenness(G)
ebc

JS <- jaccard_sim(gtc,ebc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","EB")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

#UK

data(UKfaculty,package="igraphdata")
V(UKfaculty)

gtc <- cluster_walktrap(G)
gtc

louc <- cluster_louvain(G)
louc

JS <- jaccard_sim(gtc,louc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lou")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)


labc <- cluster_label_prop(G)
labc

JS <- jaccard_sim(gtc,labc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","Lab")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

ebc <- cluster_edge_betweenness(G)
ebc

JS <- jaccard_sim(gtc,ebc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","EB")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

imc <- cluster_infomap(G)
imc

JS <- jaccard_sim(gtc,imc)
JS
print(xtable(JS, type = "latex"))
MC <- match_clusters(JS,"GT","EB")
MC
print(xtable(MC, type = "latex"))
Wmean(MC,gtc)

