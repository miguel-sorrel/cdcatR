source("nonparametric_081019.R")
source("0_Functions_CDCATcal_APM.R")

# create data
R <- 2 # 10 # Replicates
K <- 5
VAR <- 0.10 # HV in Kaplan et al. (2015)
IQ <- .70   # P(1), IQ = Low Discrimination
alpha.level <- 0.05 # 2-step LR / Wald test # 0.01, 0.05/J,...

Q5 <- q.matrix <- attributepattern(K = 5)[1:26, ][-1, ]
Q5.165 <- rbind(Q5[1:5, ], Q5[1:5, ], Q5[1:5, ], do.call("rbind", rep(list(Q5), 6))) # 05/07 Smaller item bank, no differences if it is big
TRUEmodel.165 <- c(rep(1, nrow(Q5.165)/2), rep(2, 1 + nrow(Q5.165)/2))
K <- ncol(Q5.165)
L <- 2^K
prior <- rep(1/L, L)

Ne  <- 20000 # a sufficiently large sample size for the item bank calibration
Ne2 <- 5000 # a sufficiently large sample size for the validation of the CAT
createItemBank.CDM.obj     <- createItemBank.CDM(Q5.165, IQ, VAR, TRUEmodel.165, prior, Ne, Ne2)

# compute alphas
N <- 1600*2*2
dat <- createItemBank.CDM.obj$simX_c$dat[1:N, ]
Q <- createItemBank.CDM.obj$Q
att <- createItemBank.CDM.obj$simX_c$attribute[1:N, ]
model <- createItemBank.CDM.obj$model
parm <- createItemBank.CDM.obj$item.param

TRUE.COMB <- GDINA(dat = dat, Q = Q,  
                   model = "GDINA", catprob.parm = createItemBank.CDM.obj$simX_c$catprob.parm, 
                   control = list(maxitr = 0), verbose = 0) 

GDINA.COMB <- GDINA(dat = dat, Q = Q, verbose = 0)
DINA.COMB  <- GDINA(dat = dat, Q = Q, model = "DINA", verbose = 0)
DINO.COMB  <- GDINA(dat = dat, Q = Q, model = "DINO", verbose = 0)
ACDM.COMB  <- GDINA(dat = dat, Q = Q, model = "ACDM", verbose = 0)

LR2  <-  LR_2step(GDINA.obj = GDINA.COMB) 
model.alpha <- apply(LR2$LR2.p, 1, function(x) {if (max(x, na.rm = T) > alpha.level) {
  which.max(x)}
  else {
    return(0)}})
models <- rep(0, nrow(Q))
models[which(rowSums(Q) != 1)] <- model.alpha
models.alpha <- models

model.BF <- apply(LR2$LR2.p, 1, function(x) {if (max(x, na.rm = T) > alpha.level/(sum(rowSums(Q)>1)*3)) {
  which.max(x)}
  else {
    return(0)}})
models <- rep(0, nrow(Q))
models[which(rowSums(Q) != 1)] <- model.BF
models.BF <- models

res.models.nopar <- apply(LR2$LR2.p, 1, function(x) {if (max(x, na.rm = T) > alpha.level/(sum(rowSums(Q)>1)*3)) {
  which.max(x)}
  else {
    return(sample(c(1,2), 1))}})
models <- rep(0, nrow(Q))
models[which(rowSums(Q) != 1)] <- res.models.nopar
models[models == 3] <- 2
models[models == 1] <- "AND"
models[models == 2] <- "OR"
models[models == "0"] <- sample(c("AND", "OR"), size = length(models[models == "0"]), replace = TRUE)
names(models) <- NULL
models.nopar <- models 

LR2.COMB.alpha  <- GDINA(dat = dat, Q = Q, model = models.alpha, verbose = 0)
LR2.COMB.BF  <- GDINA(dat = dat, Q = Q, model = models.BF, verbose = 0)


NPC.DINA.COMB <- NPC(data = dat, Q = Q, gate = rep("AND", nrow(Q)))
NPC.DINO.COMB <- NPC(data = dat, Q = Q, gate = rep("OR", nrow(Q)))
NPC.2LR.COMB.BF <- NPC(data = dat, Q = Q, gate = models.nopar)
gate.true <- TRUEmodel.165
gate.true[gate.true == 1] <- "AND"
gate.true[gate.true == 2] <- "OR"
NPC.TRUE <- NPC(data = dat, Q = Q, gate = gate.true)

##
rbind(
  "TRUE" = ClassRate(att, personparm(TRUE.COMB, "MAP")[, -(ncol(Q)+1)])$PCV[ncol(Q)],
  "NPC.TRUE" = ClassRate(att, NPC.TRUE$HD[, -(ncol(Q)+1)])$PCV[ncol(Q)],
  "2LR.BF"= ClassRate(att, personparm(LR2.COMB.BF, "MAP")[, -(ncol(Q)+1)])$PCV[ncol(Q)],
  "2LR.alpha" = ClassRate(att, personparm(LR2.COMB.alpha, "MAP")[, -(ncol(Q)+1)])$PCV[ncol(Q)],
  "NPC.2LR" = ClassRate(att, NPC.2LR.COMB.BF$HD[, -(ncol(Q)+1)])$PCV[ncol(Q)],
  "GDINA" = ClassRate(att, personparm(GDINA.COMB, "MAP")[, -(ncol(Q)+1)])$PCV[ncol(Q)],
"DINA" = ClassRate(att, personparm(DINA.COMB, "MAP")[, -(ncol(Q)+1)])$PCV[ncol(Q)],
"DINO" = ClassRate(att, personparm(DINO.COMB, "MAP")[, -(ncol(Q)+1)])$PCV[ncol(Q)],
"ACDM" = ClassRate(att, personparm(ACDM.COMB, "MAP")[, -(ncol(Q)+1)])$PCV[ncol(Q)],
"NPC.DINA" = ClassRate(att, NPC.DINA.COMB$HD[, -(ncol(Q)+1)])$PCV[ncol(Q)],
"NPC.DINO" = ClassRate(att, NPC.DINO.COMB$HD[, -(ncol(Q)+1)])$PCV[ncol(Q)])


