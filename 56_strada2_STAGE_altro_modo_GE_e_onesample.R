#---
#  
#output: html_document
#editor_options: 
#  chunk_output_type: console
#---

# # PACCHETTI, PARAMETRI, DATI:

# Eseguire -- librerie, WD:
# ```{r}
{
setwd(dir = "/media/beoebdue/Elements/STAGE_unimib/R_WD")
require(eegkit)
require(psych)
require(stringi)
require(igraph)
# require(RCurl)
# require(curl)
# require(edfReader)
}
# ```

# Eseguire -- parametri, IMPOSTARE COME SI DESIDERA:
# ```{r}
{
nr.fin <- 10
metodo <- "spearman"
correzione <- "holm"
TIPO.SOGLIA <- "entrambi"
CUTOFF.FISSO <- .5
SOLO.SIGNIF.TF <- TRUE
SOGLIA.SIGN.PVAL <- .05
seme.casuale <- 321
nr.repliche <- 1000
minima.clique.A = 2
massima.clique.A = 30
minima.clique.B = 2
massima.clique.B = 30
}
# ```



# # Giorno_1: ricerca cliques di B medio sui soggetti del set A 
#             ricerca cliques di A medio sui soggetti del set B, 
#             esportazione risultati
# ·
{
  
verifica.cliques.onesam <- function(nr.fin = 10, 
                               L.diL.A,
                               L.diL.B,
                               metodo = "spearman",
                               correzione = "holm",
                               TIPO.SOGLIA = "entrambi",
                               CUTOFF.FISSO = .5,
                               SOLO.SIGNIF.TF = TRUE,
                             minima.clique.A = NULL,
                             minima.clique.B = NULL,
                             massima.clique.A = NULL,
                             massima.clique.B = NULL,
                               SOGLIA.SIGN.PVAL = .05,
                               seme.casuale = 321,
                               nr.repliche = 1000){
  # # FUNZIONI:
  {
  # ```{r}
  
  # # # VECCHIE:
  # # # VECCHIE:
  # # # VECCHIE:
  
  crea.lista.vltmtx.TRIAL <- function(trial, dataset){
    dati <- dataset[dataset$trial == trial , ]
    dati <- na.omit(dati)
    soggetti <- rownames(as.matrix(which(table(dati$subject) != 0)))
    lista.vltmtx.sogg <- vector(mode = "list", length = length(soggetti))
    for (i in 1: length(soggetti)){
      soggetto.i <- soggetti[i]
      a <- dati[ dati$subject == soggetto.i , ]
      vltgmtx.a <- matrix(data = a$voltage, nrow = nlevels(x = a$channel) )
      rownames(vltgmtx.a) <- a$channel[ !duplicated(a$channel) ]
      colnames(vltgmtx.a) <- 0: 255
      acnames <- rownames(vltgmtx.a)
      idx <- c(which(acnames == "X"), which(acnames == "Y"), which(acnames == "nd"))
      vltgmtx.a <- vltgmtx.a[ -idx, ]
      lista.vltmtx.sogg[[i]] <- vltgmtx.a
      names(lista.vltmtx.sogg)[i] <- soggetto.i 
    }
    return(lista.vltmtx.sogg)
  }
  
  crea.Lfin.Lsogg <- function(Lsogg, quante.finestre){
    require(stringi)
    #Ltri.Lfin.Lsogg <- list()
    #for (i in 1: length(ListadiListe)){
    Lfin.Lsogg <- vector(mode = 'list', length = quante.finestre)
    names(Lfin.Lsogg) <- "fin_" %s+% 1: quante.finestre
    #trial.i <- ListadiListe[[i]]
    for (j in 1: length(Lsogg)){
      sogg.j <- Lsogg[[j]]
      tagli <- floor(seq(from = 0, to = dim(sogg.j)[2], length.out = quante.finestre+1))
      for (k in 1: quante.finestre){
        Fin.k <- sogg.j[, seq(from = tagli[k]+1, to = tagli[k+1])] 
        Lfin.Lsogg[[k]][[j]] <- Fin.k
        names(Lfin.Lsogg[[k]])[j] <- names(Lsogg)[j]
      }
    }
    #Ltri.Lfin.Lsogg[[i]] <- Lfin.Lsogg
    #}
    #names(Ltri.Lfin.Lsogg) <- names(ListadiListe)
    return(Lfin.Lsogg)#(Ltri.Lfin.Lsogg)
  }
  
  mediacamp.matrice.V4 <- function(Lsogg){ 
    M.mediacamp <- matrix(data = 0, 
                          nrow = nrow(Lsogg[[1]]), 
                          ncol = ncol(Lsogg[[1]])
    )
    for (j in 1: length(Lsogg)){ 
      soggetto.j.finScelta <- Lsogg[[j]] # # i=1, j=1 dim=61x64
      M.mediacamp <- M.mediacamp + soggetto.j.finScelta
    }
    M.mediacamp <- M.mediacamp / length(Lsogg)
    return(M.mediacamp)
  }
  
  cor2adj <- listaDIliste.mt.cor2adj.V4.opzno05.VERS2.BIS <- function(L.obj.corr, tipo.di.soglia, # = "fissa", "entrambi"
                                                           valore.soglia.fissa = 0.5, solo.signif, soglia.pval){
    matr.corr <- L.obj.corr$r
    matr.pvalues <- L.obj.corr$p
    distr.upper.tri <- as.vector(matr.corr[upper.tri(x = matr.corr, diag = FALSE)])
    q.rtili <- quantile(x = distr.upper.tri, probs = c(.25, .75))
    #if (tipo.di.soglia == "quartile") { soglie <- c(q.rtili[1], q.rtili[2]); print(soglie) }
    #else if (tipo.di.soglia == "fissa") { soglie <- c(-1, 1) * valore.soglia.fissa; print(soglie) }
    if (tipo.di.soglia == "fissa") { soglie <- c(-1, 1) * valore.soglia.fissa; print(soglie) }
    else { soglie <- c(min(-valore.soglia.fissa, q.rtili[1]) , max(valore.soglia.fissa, q.rtili[2])); print(soglie) }
    matr.adj <- matrix(nrow = dim(matr.corr)[1], ncol = dim(matr.corr)[2])
    for (j in 1: dim(matr.corr)[1]){
      for (k in 1: dim(matr.corr)[1]){
        if (j > k) { matr.pvalues[j, k] <- matr.pvalues[k, j] }
        if (matr.corr[j, k] < soglie[1] | matr.corr[j, k] > soglie[2]){
          if (solo.signif == TRUE){
            if (matr.pvalues[j, k] < soglia.pval) { matr.adj[j, k] <- 1 } # # attenzione: "soglia.pval" è arrivato in input alla funzione ! È una soglia POSTA PER VERIFICARE LA SIGNIFICATIVITÀ DEL p.value, non c'entra nulla col vettore: "soglie" (che valuta solo il valore di rho)
            else {matr.adj[j, k] <- 0}
          }
          else { matr.adj[j, k] <- 1 }
        }
        else { matr.adj[j, k] <- 0 }
      }
    }
    rownames(matr.adj) <- rownames(matr.corr)
    colnames(matr.adj) <- colnames(matr.corr)
    #lista.adj[[i]] <- matr.adj
    return(matr.adj)
  }
  
  ritrova <- function(lista, posizione) { return(lista[[posizione]]) }
  
  Wgg_statest <- Wgg.stattest.2smpls <- function(l1, l2){ 
    matr.medie.smpl1 <- apply(X = simplify2array(x = l1), MARGIN = 1: 2, FUN = mean)
    matr.medie.smpl2 <- apply(X = simplify2array(x = l2), MARGIN = 1: 2, FUN = mean)
    diff.matr.medie <- matr.medie.smpl1 - matr.medie.smpl2
    abs.matr.medie <- apply(X = diff.matr.medie, MARGIN = 1, FUN = abs)
    #n <- length(l1)
    #m <- length(l2)
    Wgg <- sum(abs.matr.medie[upper.tri(x = abs.matr.medie, diag = FALSE)])# * sqrt((n * m) / (n + m))
    output <- Wgg
    names(output) <- 'Wgg'
    return(output)
  }
  
  bootsXeeg <- function(set.grafi1, set.grafi2, numero.repliche){
    B <- numero.repliche # 1000		
    bb <- vector(length = B)
    c.congiunto <- append(set.grafi1, set.grafi2)
    subspl1.size <- length(set.grafi1); subspl1.size
    subspl2.size <- length(set.grafi2); subspl2.size
    for (i in 1: B){
      subspl1 <- sample(x = c.congiunto, size = subspl1.size, replace = TRUE)#; subspl1; length(subspl1)
      subspl2 <- sample(x = c.congiunto, size = subspl2.size, replace = TRUE)#; subspl2; length(subspl2)
      #names(subspl1)
      #names(subspl2)
      bb[i] <- Wgg.stattest.2smpls(l1 = subspl1, l2 = subspl2)
      if (i %in% floor(quantile(1: numero.repliche, probs = c(.25, .5, .75)))) {print(paste("rep =", floor(i)))}
    }
    Wgg.oss <- Wgg.stattest.2smpls(l1 = set.grafi1, l2 = set.grafi2)
    nonparam.pvalue <- length(bb[bb > Wgg.oss]) / length(bb); nonparam.pvalue
    names(nonparam.pvalue) <- 'nonparam.pvalue'
    return(nonparam.pvalue)
  }
  
  crea.Lfin.da.matr <- function(matrice, quante.finestre){
    require(stringi)
    Lfin <- vector(mode = 'list', length = quante.finestre)
    names(Lfin) <- "fin_" %s+% 1: quante.finestre
    #for (j in 1: length(Lsogg)){
    #  sogg.j <- Lsogg[[j]]
    tagli <- floor(seq(from = 0, to = dim(matrice)[2], length.out = quante.finestre+1))
    for (k in 1: quante.finestre){
      Fin.k <- matrice[, seq(from = tagli[k]+1, to = tagli[k+1])] 
      Lfin[[k]] <- Fin.k#[[j]] <- Fin.k
      #names(Lfin[[k]])[j] <- names(Lsogg)[j]
    }
    # }
    return(Lfin)
  }
 
  # # NUOVE:
  # # NUOVE:
  # # NUOVE:
  
  {
  # testa_una_clique <- function(clique,
  #                              graph) {
  #   # sottografo <- subgraph.edges(graph = graph, eids = clique, delete.vertices = T) # no
  #   # # i numeri assegnati ai nodi della clique corrispondono alla posizione che hanno nella relativa matrice di adiacenza
  #   {
  #     # clique
  #     # as.numeric(clique)
  #     # expand.grid(as.numeric(clique), as.numeric(clique))
  #     # colnames(adj.A)[as.numeric(clique)]
  #   }
  #   {
  #     # expand.grid(clique, clique)
  #     # as.character(clique)
  #     # class(clique)
  #     # as.list(clique)
  #   }
  #   names(clique) # estrae i nomi dei nodi della clique
  #   edges <- as.matrix(expand.grid(names(clique), names(clique))) # costruisce matrice con due colonne che indicano i nomi degli edges della clique
  #   edges2 <- edges 
  #   edges2 <- edges2[-which(edges2[, 1] == edges2[, 2]), ] # togliamo i self-edges (FP1--FP1 per esempio)
  #   
  #   sottografo <- graph_from_edgelist(el = edges2, directed = F) # questo è il grafo che seleziona gli elementi della clique, ... 
  #   # plot(graph_from_edgelist(el = edges2, directed = F))
  #   # ... ma in realtà non serve: modifichiamo la funzione che calcola la statistica test, facendogli sommare solo gli edges della clique, e vedendo se su 
  #   # questi edge c'è discriminazione
  #   return(list(matr_edges = edges2, clique = sottografo, nomi_nodi_clique = names(clique)))
  # }
  }
    
  {  # stima.pval.unbalanced <- function(set.grafi1, set.grafi2, numero.repliche){
  #   B <- numero.repliche # 1000		
  #   bb <- vector(length = B)
  #   if (class(set.grafi1) != 'list') { set.grafi1 <- list(set.grafi1) }#
  #   if (class(set.grafi2) != 'list') { set.grafi2 <- list(set.grafi2) }# entrambi son già delle liste
  #   c.congiunto <- append(set.grafi1, set.grafi2)
  #   subspl1.size <- length(set.grafi1); subspl1.size
  #   subspl2.size <- length(set.grafi2); subspl2.size
  #   # # modifica
  #   w1 <- .5 * (length(set.grafi1))^-1
  #   w2 <- .5 * (length(set.grafi2))^-1
  #   balance <- c(rep(w1, length.out = length(set.grafi1)) ,
  #                rep(w2, length.out = length(set.grafi2)))
  #   ### fine modif
  #   for (i in 1: B){
  #     subspl1 <- sample(x = c.congiunto, size = subspl1.size, replace = TRUE, prob = balance)
  #     subspl2 <- sample(x = c.congiunto, size = subspl2.size, replace = TRUE, prob = balance)
  #     {# table(names(sample(x = c.congiunto, size = subspl2.size, replace = TRUE, prob = balance)))
  #       #names(subspl1)
  #       #names(subspl2)
  #       # Wgg.stattest.2smpls(l1 = subspl1, l2 = subspl2)
  #     }
  #     bb[i] <- Wgg.stattest.2smpls(l1 = subspl1, l2 = subspl2)
  #     if (i %in% floor(quantile(1: numero.repliche, probs = c(.25, .5, .75)))) {print(paste("rep =", floor(i)))}
  #   }
  #   Wgg.oss <- Wgg.stattest.2smpls(l1 = set.grafi1, l2 = set.grafi2)
  #   nonparam.pvalue <- length(bb[bb > Wgg.oss]) / length(bb); nonparam.pvalue
  #   names(nonparam.pvalue) <- 'nonparam.pvalue'
  #   return(nonparam.pvalue)
  # }
  } # vers "stima.pval.unbalanced"
  
  stima.pval.onesample <- function(set.grafi1, set.grafi2, numero.repliche){
    B <- numero.repliche # 1000		
    bb <- vector(length = B)
    if (class(set.grafi1) != 'list') { set.grafi1 <- list(set.grafi1) }#
    if (class(set.grafi2) != 'list') { set.grafi2 <- list(set.grafi2) }# entrambi son già delle liste
    c.congiunto <- append(set.grafi1, set.grafi2)
    subspl1.size <- length(set.grafi1); subspl1.size
    subspl2.size <- length(set.grafi2); subspl2.size
    # # modifica
    w1 <- .5 * (length(set.grafi1))^-1
    w2 <- .5 * (length(set.grafi2))^-1
    balance <- c(rep(w1, length.out = length(set.grafi1)) ,
                 rep(w2, length.out = length(set.grafi2)))
    ### fine modif
    for (i in 1: B){
      # subspl1 <- sample(x = c.congiunto, size = subspl1.size, replace = TRUE, prob = balance)
      subspl2 <- sample(x = c.congiunto, size = subspl2.size, replace = TRUE, prob = balance)# solo uno
      {# table(names(sample(x = c.congiunto, size = subspl2.size, replace = TRUE, prob = balance)))
        #names(subspl1)
        #names(subspl2)
        # Wgg.stattest.2smpls(l1 = subspl1, l2 = subspl2)
      }
      bb[i] <- Wgg.stattest.2smpls(l1 = set.grafi1, # # # # # # # # # # di qua gli passiamo il set che aveva SOLAMENTE LA MATRICE MEDIA ... HA SENSO? 
                                   l2 = subspl2)
      if (i %in% floor(quantile(1: numero.repliche, probs = c(.25, .5, .75)))) {print(paste("rep =", floor(i)))}
    }
    Wgg.oss <- Wgg.stattest.2smpls(l1 = set.grafi1, l2 = set.grafi2)
    nonparam.pvalue <- length(bb[bb >= Wgg.oss]) / length(bb); nonparam.pvalue
    names(nonparam.pvalue) <- 'nonparam.pvalue'
    return(nonparam.pvalue)
  }


  
  }
  
  # # APPLICAZIONE PROCESSO:
  # ·
  if (class(L.diL.A) != 'list') { L.diL.A <- list(L.diL.A)}
  if (class(L.diL.B) != 'list') { L.diL.B <- list(L.diL.B)}
  {
  #   1) Medie:
  # ```{r}
  # L.diMedie.A <- list()
  # for (i in 1: length(L.diL.A)){
  #   cuTRA <- L.diL.A[i]
  #   Media.A <- apply(X = simplify2array(x = cuTRA), MARGIN = 1: 2, FUN = mean)
  #   L.diMedie.A[[i]] <- Media.A
  # }
  }
  Media.TaskA <- apply(X = simplify2array(x = L.diL.A
  ), MARGIN = 1: 2, FUN = mean)
  # ```
  
  # ```{r}
  {
  # L.diMedie.B <- list()
  # for (i in 1: length(L.diL.B)){
  #   cuTRB <- L.diL.B[i]
  #   Media.B <- apply(X = simplify2array(x = cuTRB), MARGIN = 1: 2, FUN = mean)
  #   L.diMedie.B[[i]] <- Media.B
  # }
  }
  Media.TaskB <- apply(X = simplify2array(x = L.diL.B
  ), MARGIN = 1: 2, FUN = mean)
  # ```
  
  #   2) traspongo: (NON taglio finestre)
  # ```{r}
  {
  # Lfin.MTaskA <- crea.Lfin.da.matr(matrice = Media.TaskA, quante.finestre = nr.fin)
  # Lfin.MTaskB <- crea.Lfin.da.matr(matrice = Media.TaskB, quante.finestre = nr.fin)
  # # ```
  # # ```{r}
  # TLfin.MTaskA <- lapply(X = Lfin.MTaskA, FUN = t)
  # TLfin.MTaskB <- lapply(X = Lfin.MTaskB, FUN = t)
  }
  T.MTaskA <- t(Media.TaskA)#lapply(X = Media.TaskA, FUN = t)
  T.MTaskB <- t(Media.TaskB)#lapply(X = Media.TaskB, FUN = t)
  # ```

  #   3) Da matrici delle medie campionarie a CORR (con correzione test multipli)
  # ```{r}
  {
  # Lfin.corr.A <- lapply(X = TLfin.MTaskA, FUN =  
  #                           corr.test, method = metodo, adjust = correzione, ci = FALSE)
  # Lfin.corr.B <- lapply(X = TLfin.MTaskB, FUN = 
  #                           corr.test, method = metodo, adjust = correzione, ci = FALSE)
  }
  corr.A <- corr.test(T.MTaskA, method = metodo, adjust = correzione, ci = FALSE)
  corr.B <- corr.test(T.MTaskB, method = metodo, adjust = correzione, ci = FALSE)
  # ```

  #   4) da CORR a ADJ 
  # ```{r}
  {
  # Lfin.adj.A <- lapply(X = Lfin.corr.A, FUN = listaDIliste.mt.cor2adj.V4.opzno05.VERS2.BIS,
  #                        tipo.di.soglia = TIPO.SOGLIA, valore.soglia.fissa = CUTOFF.FISSO,
  #                        solo.signif = SOLO.SIGNIF.TF, soglia.pval = SOGLIA.SIGN.PVAL)
  # Lfin.adj.B <- lapply(X = Lfin.corr.B, FUN = listaDIliste.mt.cor2adj.V4.opzno05.VERS2.BIS, 
  #                        tipo.di.soglia = TIPO.SOGLIA, valore.soglia.fissa = CUTOFF.FISSO,
  #                        solo.signif = SOLO.SIGNIF.TF, soglia.pval = SOGLIA.SIGN.PVAL)
  }
  identical(x = cor2adj, y = listaDIliste.mt.cor2adj.V4.opzno05.VERS2.BIS) # solo un memorandum
  adj.A <- cor2adj(L.obj.corr = corr.A, 
                   tipo.di.soglia = TIPO.SOGLIA, 
                   valore.soglia.fissa = CUTOFF.FISSO, 
                   solo.signif = SOLO.SIGNIF.TF, 
                   soglia.pval = SOGLIA.SIGN.PVAL)
  adj.B <- cor2adj(L.obj.corr = corr.B, 
                   tipo.di.soglia = TIPO.SOGLIA, 
                   valore.soglia.fissa = CUTOFF.FISSO, 
                   solo.signif = SOLO.SIGNIF.TF, 
                   soglia.pval = SOGLIA.SIGN.PVAL)
  adj.A.diag0 <- adj.A
  diag(adj.A.diag0) <- 0
  adj.B.diag0 <- adj.B
  diag(adj.B.diag0) <- 0
  
  # ```
######## AGGIUNTA:
  # # # # # # # # AGGIUNTA:
  # # # # # # # # AGGIUNTA:
  # # # # # # # # AGGIUNTA:

  # convertire in igraph
  igr.adj.A <- graph.adjacency(adjmatrix = adj.A.diag0, mode = "undirected", weighted = NULL, diag = F)
  igr.adj.B <- graph.adjacency(adjmatrix = adj.B.diag0, mode = "undirected", weighted = NULL, diag = F)

  # sistemazioni varie
  if (is.null(massima.clique.A)) {
    massima.clique.A <- clique_num(igr.adj.A) 
    warning("your 'massima.clique.A' parameter is NULL; it has been changed to " %s+% clique_num(igr.adj.A) %s+% ", that is the size of the largest clique(s) for group A")
  } else if (clique_num(igr.adj.A) < massima.clique.A) {
    massima.clique.A <- clique_num(igr.adj.A) 
    warning("the size of the largest clique(s) for group A is " %s+% clique_num(igr.adj.A) %s+% "; your 'massima.clique.A' parameter is higher than this number, it has been changed to " %s+% clique_num(igr.adj.A))
  }
  if (is.null(massima.clique.B)) {
    massima.clique.B <- clique_num(igr.adj.B) 
    warning("your 'massima.clique.B' parameter is NULL; it has been changed to " %s+% clique_num(igr.adj.B) %s+% ", that is the size of the largest clique(s) for group B")
  } else if (clique_num(igr.adj.B) < massima.clique.B) {
    massima.clique.B <- clique_num(igr.adj.B) 
    warning("the size of the largest clique(s) for group B is " %s+% clique_num(igr.adj.B) %s+% "; your 'massima.clique.B' parameter is higher than this number, it has been changed to " %s+% clique_num(igr.adj.B))
  }
  
  if (is.null(minima.clique.A)) {
    minima.clique.A <- massima.clique.A
    warning("argument not provided to 'minima.clique.A'; it has been set to " %s+% massima.clique.A %s+% ", that is equal to 'massima.clique.A'")
  }
  if (is.null(minima.clique.B)) {
    minima.clique.B <- massima.clique.B
    warning("argument not provided to 'minima.clique.B'; it has been set to " %s+% massima.clique.B %s+% ", that is equal to 'massima.clique.B'")
  }
  
  # individuare cliques massime di A e di B
  cliques.A <- max_cliques(graph = igr.adj.A, 
                            min = minima.clique.A,
                            max = massima.clique.A)
  cliques.B <- max_cliques(graph = igr.adj.B, 
                            min = minima.clique.B,
                            max = massima.clique.B)
  
  # ottenere i nomi di riga / colonna che costituiscono le cliques di A
  nomi.clq.A <- simplify2array(lapply(X = seq_along(cliques.A),
                                      FUN = function(i, Lclq) { names(Lclq[[i]]) },
                                      Lclq = cliques.A))
  # ottenere i nomi di riga / colonna che costituiscono le cliques di B
  nomi.clq.B <- simplify2array(lapply(X = seq_along(cliques.B),
                                      FUN = function(i, Lclq) { names(Lclq[[i]]) },
                                      Lclq = cliques.B))
  
  # # 1. Quanto A somiglia a B sugli archi in cui B ha cliques ?
  # adj.B[nomi.clq.B[, 1], nomi.clq.B[, 1]]
  # 
  # cose da fare:
  #   i· fare un grafo per ogni soggetto di A (corr -> adj) e selezionarne gli edge con nomi.clq.B
  #   ii· fare la valutazione di dissomiglianza per ogni clique di B
  #   iii· esportarne i risultati 
  # e stessa cosa per B con A
  
  # # (i) un grafo per ogni soggetto di A ...
  # L.diL.A; length(L.diL.A)
  T.L.diL.A <- lapply(L.diL.A, t)
  Lcorr.sA <- lapply(X = T.L.diL.A, 
                     FUN = corr.test, method = metodo, adjust = correzione, ci = FALSE)
  Ladj.sA <- lapply(X = Lcorr.sA, 
                    FUN = cor2adj,
                    tipo.di.soglia = TIPO.SOGLIA, 
                    valore.soglia.fissa = CUTOFF.FISSO, 
                    solo.signif = SOLO.SIGNIF.TF, 
                    soglia.pval = SOGLIA.SIGN.PVAL)
  # ... e di B
  T.L.diL.B <- lapply(L.diL.B, t)
  Lcorr.sB <- lapply(X = T.L.diL.B, 
                     FUN = corr.test, method = metodo, adjust = correzione, ci = FALSE)
  Ladj.sB <- lapply(X = Lcorr.sB, 
                    FUN = cor2adj,
                    tipo.di.soglia = TIPO.SOGLIA, 
                    valore.soglia.fissa = CUTOFF.FISSO, 
                    solo.signif = SOLO.SIGNIF.TF, 
                    soglia.pval = SOGLIA.SIGN.PVAL)
  
  #print(class(nomi.clq.B)); print(class(nomi.clq.A))
  #print(nomi.clq.B)
  #print(nomi.clq.A)
  len.B <- length(nomi.clq.B)#dim(nomi.clq.B)[2]
  len.A <- length(nomi.clq.A)#dim(nomi.clq.A)[2]
  
  # # (ii) valutazione dissomiglianza per ogni soggetto di A da ciascuna clique di B ...
  dissim.AfromB <- lapply(1: len.B,
                          FUN = function(i, Ls, ADJ.to.refer, nomi.clq, nr.repliche) { 
                            clique <- nomi.clq[[i]]
                            print('clique_' %s+% i %s+% '_A-B')
                            Ls.select <- lapply(Ls, function(s, cl) { return(s[cl, cl]) },
                                                cl = clique)
                            # diagonale =1 non importa, tanto la funzione Wgg.statest ignora la diagonale
                            ADJ.ref.select <- ADJ.to.refer[clique , clique]
                            Wgg.oss <- Wgg_statest(l1 = ADJ.ref.select, l2 = Ls.select) # stat.test.oss
                            # pvalue stimato
                            stima.p.value <- stima.pval.onesample(set.grafi1 = ADJ.ref.select, 
                                                                   set.grafi2 = Ls.select, 
                                                                   numero.repliche = nr.repliche)
                            stringa.nomi <- paste0(clique, collapse = ',')
                            return(data.frame(clq = stringa.nomi,
                                              clq_size = length(clique),
                                              Wgg_oss = Wgg.oss,
                                              estim_pvalue = stima.p.value, row.names = NULL))
                          }, 
                          Ls = Ladj.sA,
                          ADJ.to.refer = adj.B.diag0,
                          nomi.clq = nomi.clq.B,
                          nr.repliche = nr.repliche)
  names(dissim.AfromB) <- paste0("setAcliqB__Clq", 1: length(dissim.AfromB))
  # # (ii) ... e per ogni soggetto di B da ciascuna clique di A
  dissim.BfromA <- lapply(1: len.A,
                          FUN = function(i, Ls, ADJ.to.refer, nomi.clq, nr.repliche) { 
                            clique <- nomi.clq[[i]]
                            print('clique_' %s+% i %s+% '_B-A')
                            Ls.select <- lapply(Ls, function(s, cl) { return(s[cl, cl]) },
                                                cl = clique)
                            # diagonale =1 non importa, tanto la funzione Wgg.statest ignora la diagonale
                            ADJ.ref.select <- ADJ.to.refer[clique , clique]
                            Wgg.oss <- Wgg_statest(l1 = ADJ.ref.select, l2 = Ls.select) # stat.test.oss
                            # pvalue stimato
                            stima.p.value <- stima.pval.onesample(set.grafi1 = ADJ.ref.select, 
                                                                   set.grafi2 = Ls.select, 
                                                                   numero.repliche = nr.repliche)
                            stringa.nomi <- paste0(clique, collapse = ',')
                            return(data.frame(clq = stringa.nomi,
                                              clq_size = length(clique),
                                              Wgg_oss = Wgg.oss,
                                              estim_pvalue = stima.p.value, row.names = NULL))
                          }, 
                          Ls = Ladj.sB,
                          ADJ.to.refer = adj.A.diag0,
                          nomi.clq = nomi.clq.A,
                          nr.repliche = nr.repliche)
  names(dissim.BfromA) <- paste0("setBcliqA__Clq", 1: length(dissim.BfromA))
  
  print(dissim.AfromB)
  print(dissim.BfromA)
  # # (iii) esportare
  return(append(dissim.AfromB, dissim.BfromA))
  
  # # # # # # # # # #
  # # # # # # # # # # FINE
  # # # # # # # # # # 
######### FINE

  # ```
}
# ·
# # 
# ·
{
# per lanciare questo codice:
#
# eseguire il punto (0) e il punto (1) del codice 28_STAGE_serv, 
#
# OPPURE importare gli oggetti: "DATASET_2_cuTR_03_07_11.RData" e "DATASET_2_cuTR_04_08_12.RData":
# load("DATASET_2_cuTR_03_07_11.RData") # togli commento iniziale da questa riga
# load("DATASET_2_cuTR_04_08_12.RData") # togli commento iniziale da questa riga
# e lanciare:
# L.di.L.1.1 <- list(cuTR03, cuTR07, cuTR11) # togli commento iniziale da questa riga
# L.di.L.2.1 <- list(cuTR04, cuTR08, cuTR12) # togli commento iniziale da questa riga
#
# OPPURE passare alla funzione direttamente le matrici-medie-elementwise contenute 
# negli oggetti: "DATASET_2_Media_03_07_11.RData" e "DATASET_2_Media_04_08_12.RData"
# (dentro la funzione farà la media elementwise dell'unica matrice )
 
# load("DATASET_2_Media_03_07_11.RData")
# load("DATASET_2_Media_04_08_12.RData")
} # <~~ no! devi caricare oggetto r con LsoggA e LsoggC

L.di.L.1.1 <- Ltrials.Lsogg.A$`111`# Lsogg.A # 
L.di.L.2.1 <- Ltrials.Lsogg.C$`111`# Lsogg.C # è solo trial 15

VKLK.onesam <- verifica.cliques.onesam(nr.fin = nr.fin, 
                          L.diL.A = L.di.L.1.1,
                          L.diL.B = L.di.L.2.1, 
                          metodo = metodo,
                          correzione = correzione,
                          TIPO.SOGLIA = TIPO.SOGLIA,
                          CUTOFF.FISSO = CUTOFF.FISSO,
                          SOLO.SIGNIF.TF = SOLO.SIGNIF.TF,
                          SOGLIA.SIGN.PVAL = SOGLIA.SIGN.PVAL,
                          seme.casuale = seme.casuale,
                          nr.repliche = nr.repliche,
                     minima.clique.A = minima.clique.A,
                     massima.clique.A = massima.clique.A,
                     minima.clique.B = minima.clique.B,
                     massima.clique.B = massima.clique.B)

VKLK.onesam
# ·
# # 

# # # # # # # # 

# # # # # # # # STAMPA su file:
{
  # tabella <- data.frame(statest_lp, pvalue_lp, statest_betw, pvalue_betw, 
  #                       statest_clust, pvalue_clust, statest_modul, pvalue_modul)
  data.frame(VKLK.onesam)
  tabella <- as.data.frame(t(data.frame(lapply(X = VKLK.onesam, FUN = t))))
  class(tabella[, 1])
  class(tabella[, 2])
  class(tabella[, 3])
  class(tabella[, 4])
  # class(tabella[, 1])
  tabella[, 2] <- as.numeric(as.character(tabella[, 2]))#[f]as.numeric(tabella[, 2])
  tabella[, 3] <- as.numeric(as.character(tabella[, 3]))#as.numeric(tabella[, 3])
  tabella[, 4] <- as.numeric(as.character(tabella[, 4]))#as.numeric(tabella[, 4])
  
  require(lubridate)
  tempo <- stri_replace_all(str = lubridate::now(), regex = ':', replacement = '.')
  tempo <- stri_replace_all(str = tempo, regex = ' ', replacement = '_')
  file.create("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/" %s+% 
                tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_risult_cliques.csv")
  write.table(x = tabella, 
              file = "/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/" %s+% 
                        tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_risult_cliques.csv", 
              append = F, 
              quote = T, sep = ';', dec = '.',
              row.names = T)#as.character(prMinTrials))#F)
}

#   7) sistemare a mano l'intestazione del file appena creato [aggiungere uno spazio all'inizio]
# 

#   8) Registra parametri:
#
{
  zz <- file("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/" %s+% 
                tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_risult_cliques.csv",
             "a")  # open an output file connection
  writeLines(text = ("Parametri: nr.finestre = " %s+% nr.fin %s+% 
                       ", correl: " %s+% metodo %s+% 
                       ", correzione: " %s+% correzione %s+% 
                       ", singif.correl = " %s+% SOGLIA.SIGN.PVAL %s+% 
                       ", seme.casuale = " %s+% seme.casuale %s+% 
                       ",  nr.repliche = " %s+% nr.repliche  %s+% 
                     ", minima.clique.A = "  %s+% minima.clique.A  %s+%  
                     ", massima.clique.A = "  %s+% massima.clique.A  %s+%   
                     ", minima.clique.B = "  %s+% minima.clique.B  %s+%   
                     ", massima.clique.B = "  %s+% massima.clique.B),
             con = zz)
  close(zz)
}


# # # # # PER TUTTI I TRIALS:
# # # # # PER TUTTI I TRIALS: calcola dissimilarità
# # # # # PER TUTTI I TRIALS:
#  da codice 21:
{
  crea.lista.vltmtx.TRIAL <- function(trial, dataset){
    dati <- dataset[dataset$trial == trial , ]
    dati <- na.omit(dati)
    soggetti <- rownames(as.matrix(which(table(dati$subject) != 0)))
    lista.vltmtx.sogg <- vector(mode = "list", length = length(soggetti))
    for (i in 1: length(soggetti)){
      soggetto.i <- soggetti[i]
      a <- dati[ dati$subject == soggetto.i , ]
      vltgmtx.a <- matrix(data = a$voltage, nrow = nlevels(x = a$channel) )
      rownames(vltgmtx.a) <- a$channel[ !duplicated(a$channel) ]
      colnames(vltgmtx.a) <- 0: 255
      acnames <- rownames(vltgmtx.a)
      idx <- c(which(acnames == "X"), which(acnames == "Y"), which(acnames == "nd"))
      vltgmtx.a <- vltgmtx.a[ -idx, ]
      lista.vltmtx.sogg[[i]] <- vltgmtx.a
      names(lista.vltmtx.sogg)[i] <- soggetto.i 
    }
    return(lista.vltmtx.sogg)
  }
  
  load(file = "wrkspc_DATI_eegV2A.RData") 
  load(file = "wrkspc_DATI_eegV2C.RData") 
  
  # # le trials sono:
  table(eeg.V2.A$trial[eeg.V2.A$subject == 'co3a0000448'])
  prMinTrials <- cbind(15,  19,  29,  35,  45,  49,  53,  63,  67,  83,  87,  95, 103, 111)
  
  # # lista di trials di soggetti, gruppo A e C
  Ltrials.Lsogg.A <- apply(X = prMinTrials, MARGIN = 2, FUN = crea.lista.vltmtx.TRIAL, dataset = eeg.V2.A)
  Ltrials.Lsogg.C <- apply(X = prMinTrials, MARGIN = 2, FUN = crea.lista.vltmtx.TRIAL, dataset = eeg.V2.C)
  names(Ltrials.Lsogg.A) <- as.character(prMinTrials)
  Ltrials.Lsogg.A
  names(Ltrials.Lsogg.C) <- as.character(prMinTrials)
  Ltrials.Lsogg.C
}
Ltrials.VKLK.onesam <- lapply(X = 1: length(prMinTrials), 
                       FUN = function(i, 
                                      Ltr1,
                                      Ltr2,
                                      nr.fin, 
                                      metodo,
                                      correzione,
                                      TIPO.SOGLIA ,
                                      CUTOFF.FISSO,
                                      SOLO.SIGNIF.TF,
                                      SOGLIA.SIGN.PVAL,
                                      seme.casuale,
                                      nr.repliche,
                                      minima.clique.A,
                                      massima.clique.A,
                                      minima.clique.B,
                                      massima.clique.B) { 
                         TR1i <- Ltr1[[i]]
                         TR2i <- Ltr2[[i]]
                         print('#==============================================# #==============================================#')
                         print('#==============================================# #==============================================#')
                         print("TRIAL:: " %s+% names(Ltr1)[i])
                         VKLK.onesam_i <- verifica.cliques.onesam(nr.fin = nr.fin, 
                                                    L.diL.A = TR1i,
                                                    L.diL.B = TR2i, 
                                                    metodo = metodo,
                                                    correzione = correzione,
                                                    TIPO.SOGLIA = TIPO.SOGLIA,
                                                    CUTOFF.FISSO = CUTOFF.FISSO,
                                                    SOLO.SIGNIF.TF = SOLO.SIGNIF.TF,
                                                    SOGLIA.SIGN.PVAL = SOGLIA.SIGN.PVAL,
                                                    seme.casuale = seme.casuale,
                                                    nr.repliche = nr.repliche,
                                                    minima.clique.A = minima.clique.A,
                                                    massima.clique.A = massima.clique.A,
                                                    minima.clique.B = minima.clique.B,
                                                    massima.clique.B = massima.clique.B)
                         
                         },
                       Ltr1 = Ltrials.Lsogg.A,
                       Ltr2 = Ltrials.Lsogg.C,
                       nr.fin = nr.fin, 
                       metodo = metodo,
                       correzione = correzione,
                       TIPO.SOGLIA = TIPO.SOGLIA,
                       CUTOFF.FISSO = CUTOFF.FISSO,
                       SOLO.SIGNIF.TF = SOLO.SIGNIF.TF,
                       SOGLIA.SIGN.PVAL = SOGLIA.SIGN.PVAL,
                       seme.casuale = seme.casuale,
                       nr.repliche = nr.repliche,
                       minima.clique.A = minima.clique.A,
                       massima.clique.A = massima.clique.A,
                       minima.clique.B = minima.clique.B,
                       massima.clique.B = massima.clique.B)
names(Ltrials.VKLK.onesam) <- paste0("Trial_", prMinTrials)


Ltrials.VKLK.onesam$Trial_111








# scrivi risult di tutte le trials
# scrivi risult di tutte le trials
# scrivi risult di tutte le trials

# require(parallel)
# # Use the detectCores() function to find the number of cores in system
# no_cores <- detectCores()
# # Setup cluster
# clust <- makeCluster(no_cores) #This line will take time

lapply(   X = 1: dim(prMinTrials)[2], 
       FUN = function(i, 
                      Ltr.VKLK.onesam, 
                      prMinTr, 
                      nr.repliche,
                      nr.fin
                      , metodo
                      , correzione
                      , SOGLIA.SIGN.PVAL
                      , seme.casuale 
                      , minima.clique.A
                      , massima.clique.A  
                      , minima.clique.B  
                      , massima.clique.B ) {
         require(stringi)
         VKLK.onesam <- Ltr.VKLK.onesam[[i]]
         trial <- prMinTr[, i]
         {

           # tabella <- data.frame(statest_lp, pvalue_lp, statest_betw, pvalue_betw, 
           #                       statest_clust, pvalue_clust, statest_modul, pvalue_modul)
           data.frame(VKLK.onesam)
           tabella <- as.data.frame(t(data.frame(lapply(X = VKLK.onesam, FUN = t))))
           class(tabella[, 1])
           class(tabella[, 2])
           class(tabella[, 3])
           class(tabella[, 4])
           # class(tabella[, 1])
           tabella[, 2] <- as.numeric(as.character(tabella[, 2]))#[f]as.numeric(tabella[, 2])
           tabella[, 3] <- as.numeric(as.character(tabella[, 3]))#as.numeric(tabella[, 3])
           tabella[, 4] <- as.numeric(as.character(tabella[, 4]))#as.numeric(tabella[, 4])
           
           require(lubridate)
           tempo <- stri_replace_all(str = lubridate::now(), regex = ':', replacement = '.')
           tempo <- stri_replace_all(str = tempo, regex = ' ', replacement = '_')
           file.create("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
                         tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_Trial"  %s+%  trial %s+% "_risult_cliques.csv")
           write.table(x = tabella, 
                       file = "/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
                         tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_Trial"  %s+%  trial %s+% "_risult_cliques.csv", 
                       append = F, 
                       quote = T, sep = ';', dec = '.',
                       row.names = T)#as.character(prMinTrials))#F)
         }
         #   7) sistemare a mano l'intestazione del file appena creato [aggiungere uno spazio all'inizio]
         # 
         #   8) Registra parametri:
         #
         {
           zz <- file("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
                        tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_Trial"  %s+%  trial %s+% "_risult_cliques.csv",
                      "a")  # open an output file connection
           writeLines(text = ("Parametri: nr.finestre = " %s+% nr.fin %s+% 
                                ", correl: " %s+% metodo %s+% 
                                ", correzione: " %s+% correzione %s+% 
                                ", singif.correl = " %s+% SOGLIA.SIGN.PVAL %s+% 
                                ", seme.casuale = " %s+% seme.casuale %s+% 
                                ", nr.repliche = " %s+% nr.repliche  %s+% 
                                ", minima.clique.A = "  %s+% minima.clique.A  %s+%  
                                ", massima.clique.A = "  %s+% massima.clique.A  %s+%   
                                ", minima.clique.B = "  %s+% minima.clique.B  %s+%   
                                ", massima.clique.B = "  %s+% massima.clique.B),
                      con = zz)
           close(zz)
         }
       },
       Ltr.VKLK.onesam = Ltrials.VKLK.onesam,
       prMinTr = prMinTrials,
       nr.repliche = nr.repliche,
       nr.fin = nr.fin
         , metodo =  metodo 
         , correzione =  correzione 
         , SOGLIA.SIGN.PVAL =  SOGLIA.SIGN.PVAL 
         , seme.casuale =  seme.casuale 
         , minima.clique.A =   minima.clique.A   
         , massima.clique.A =   massima.clique.A    
         , minima.clique.B =   minima.clique.B    
         , massima.clique.B =   massima.clique.B)

# stopCluster(clust)


# # # # REMEMBER: 
# we passed Alcoholic subjects as first list, so in the exported file "A" group is "Alcoholic"
# as second list control subjects, so in the exported file "B" label means "Control"
# we should change the files

}



# # # # # # # # # # # # # # # ESPLORAZIONE RISULTATI:
# # # # # # # # # # # # # # # ESPLORAZIONE RISULTATI:

# crea unico dataframe
{
listabella_onesam <-   lapply(   X = 1: dim(prMinTrials)[2], 
          FUN = function(i, 
                         Ltr.VKLK.onesam, 
                         prMinTr, 
                         nr.repliche,
                         nr.fin
                         , metodo
                         , correzione
                         , SOGLIA.SIGN.PVAL
                         , seme.casuale 
                         , minima.clique.A
                         , massima.clique.A  
                         , minima.clique.B  
                         , massima.clique.B ) {
            require(stringi)
            VKLK.onesam <- Ltr.VKLK.onesam[[i]]
            trial <- prMinTr[, i]
              # tabella <- data.frame(statest_lp, pvalue_lp, statest_betw, pvalue_betw, 
              #                       statest_clust, pvalue_clust, statest_modul, pvalue_modul)
              data.frame(VKLK.onesam)
              tabella <- as.data.frame(t(data.frame(lapply(X = VKLK.onesam, FUN = t))))
              lis <- stri_split_regex(str = rownames(tabella), pattern = '__', simplify = T)
              col1 <- stri_replace(str = lis[, 1], fixed = "cliq", replacement = "_")
              col1 <- stri_replace(str = col1, fixed = "set", replacement = "")
              tabella <- data.frame(tabella, da_a_ = col1, quale_clique = lis[, 2], row.names = "trial_"  %s+%  trial  %s+% '_'  %s+%  row.names(tabella))
              {              # class(tabella[, 1])
              # class(tabella[, 2])
              # class(tabella[, 3])
              # class(tabella[, 4])
              # class(tabella[, 1])
              }
              tabella[, 2] <- as.numeric(as.character(tabella[, 2]))#[f]as.numeric(tabella[, 2])
              tabella[, 3] <- as.numeric(as.character(tabella[, 3]))#as.numeric(tabella[, 3])
              tabella[, 4] <- as.numeric(as.character(tabella[, 4]))#as.numeric(tabella[, 4])
              {
            #   require(lubridate)
            #   tempo <- stri_replace_all(str = lubridate::now(), regex = ':', replacement = '.')
            #   tempo <- stri_replace_all(str = tempo, regex = ' ', replacement = '_')
            #   file.create("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
            #                 tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_Trial"  %s+%  trial %s+% "_risult_cliques.csv")
            #   write.table(x = tabella, 
            #               file = "/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
            #                 tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_Trial"  %s+%  trial %s+% "_risult_cliques.csv", 
            #               append = F, 
            #               quote = T, sep = ';', dec = '.',
            #               row.names = T)#as.character(prMinTrials))#F)
            # }
            # #   7) sistemare a mano l'intestazione del file appena creato [aggiungere uno spazio all'inizio]
            # # 
            # #   8) Registra parametri:
            # #
            # {
            #   zz <- file("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
            #                tempo %s+% "_n.rep" %s+% nr.repliche %s+% "_Trial"  %s+%  trial %s+% "_risult_cliques.csv",
            #              "a")  # open an output file connection
            #   writeLines(text = ("Parametri: nr.finestre = " %s+% nr.fin %s+% 
            #                        ", correl: " %s+% metodo %s+% 
            #                        ", correzione: " %s+% correzione %s+% 
            #                        ", singif.correl = " %s+% SOGLIA.SIGN.PVAL %s+% 
            #                        ", seme.casuale = " %s+% seme.casuale %s+% 
            #                        ", nr.repliche = " %s+% nr.repliche  %s+% 
            #                        ", minima.clique.A = "  %s+% minima.clique.A  %s+%  
            #                        ", massima.clique.A = "  %s+% massima.clique.A  %s+%   
            #                        ", minima.clique.B = "  %s+% minima.clique.B  %s+%   
            #                        ", massima.clique.B = "  %s+% massima.clique.B),
            #              con = zz)
            #   close(zz)
            # }
              }
              return(tabella)
              
          },
          Ltr.VKLK.onesam = Ltrials.VKLK.onesam,
          prMinTr = prMinTrials,
          nr.repliche = nr.repliche,
          nr.fin = nr.fin
          , metodo =  metodo 
          , correzione =  correzione 
          , SOGLIA.SIGN.PVAL =  SOGLIA.SIGN.PVAL 
          , seme.casuale =  seme.casuale 
          , minima.clique.A =   minima.clique.A   
          , massima.clique.A =   massima.clique.A    
          , minima.clique.B =   minima.clique.B    
          , massima.clique.B =   massima.clique.B )
}  
RESULT_onesam <- do.call(rbind, listabella_onesam) 
# già che ci siamo esportiamolo
{
  require(lubridate)
  tempo <- stri_replace_all(str = lubridate::now(), regex = ':', replacement = '.')
  tempo <- stri_replace_all(str = tempo, regex = ' ', replacement = '_')
  file.create("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
                "TABELLONE_RESULT_onesam_" %s+% tempo %s+% "_n.rep" %s+% nr.repliche %s+% ".csv")
  write.table(x = RESULT_onesam,
              file = "/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
                "TABELLONE_RESULT_onesam_" %s+% tempo %s+% "_n.rep" %s+% nr.repliche %s+% ".csv",
              append = F,
              quote = T, sep = ';', dec = '.',
              row.names = T)#as.character(prMinTrials))#F)
}
{
  zz <- file("/media/beoebdue/Elements/STAGE_unimib/DA MANDARE PER CONCLUDERE expansion/trials_onesam/" %s+% 
               "TABELLONE_RESULT_onesam_" %s+% tempo %s+% "_n.rep" %s+% nr.repliche %s+% ".csv",
             "a")  # open an output file connection
  writeLines(text = ("Parametri: nr.finestre = " %s+% nr.fin %s+%
                       ", correl: " %s+% metodo %s+%
                       ", correzione: " %s+% correzione %s+%
                       ", singif.correl = " %s+% SOGLIA.SIGN.PVAL %s+%
                       ", seme.casuale = " %s+% seme.casuale %s+%
                       ", nr.repliche = " %s+% nr.repliche  %s+%
                       ", minima.clique.A = "  %s+% minima.clique.A  %s+%
                       ", massima.clique.A = "  %s+% massima.clique.A  %s+%
                       ", minima.clique.B = "  %s+% minima.clique.B  %s+%
                       ", massima.clique.B = "  %s+% massima.clique.B),
             con = zz)
  close(zz)
}
  
Ltrials.VKLK.onesam$Trial_15$setAcliqB__Clq1$clq

table(RESULT_onesam$clq)
density(plot(table(RESULT_onesam$clq)))
table(RESULT_onesam$clq)[order(table(RESULT_onesam$clq), decreasing = F)]


RESULT_onesam.ge.15 <- RESULT_onesam[which(RESULT_onesam$clq_size >= 15), ]
table(RESULT_onesam.ge.15$clq)[table(RESULT_onesam.ge.15$clq) != 0]



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # TABELLA DI CONTINGENZA PVALUES A_B B_A E ISTOGRAMMI

table(ifelse(RESULT_onesam$estim_pvalue <= .05, 'signif', 'non_signif'), RESULT_onesam$da_a_)
tabellina.signif <- table(ifelse(RESULT_onesam$estim_pvalue <= .05, 'signif', 'non_signif'), RESULT_onesam$da_a_)
chisq.test(tabellina.signif)


par(mfrow = c(1, 2))
hist(RESULT_onesam$estim_pvalue[RESULT_onesam$da_a_ == 'A_B'], breaks = 30); abline(v = .05, col = 'red')
hist(RESULT_onesam$estim_pvalue[RESULT_onesam$da_a_ == 'B_A'], breaks = 30); abline(v = .05, col = 'red')
hist(RESULT_onesam$estim_pvalue[RESULT_onesam$da_a_ == 'A_B'], breaks = 90, main = "A_B p-values", xlab = ''); abline(v = .05, col = 'red')#; text(x = .05,  pos = 2, '0.05')#points(x = c(.05, 0), add = T, pch = 'x', col = 'red')
hist(RESULT_onesam$estim_pvalue[RESULT_onesam$da_a_ == 'B_A'], breaks = 90, main = "B_A p-values", xlab = '', ylim = c(0, 400)); abline(v = .05, col = 'red')
par(mfrow = c(1, 1))

{1-(46/(46+488))
329/(230+329)
488+46
488/534}








