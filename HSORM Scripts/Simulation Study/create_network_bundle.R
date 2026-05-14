# Minimal deps: igraph (local only), Matrix
library(igraph)
library(dplyr)
library(Matrix)

# INPUTS you already have locally:
# - stepwise_network_list : list of igraph graphs by step (your code creates this)
# - npi_trial_data        : data.frame with columns c("npi","trial_step","HospitalID")

outdir<-"/drives/drive1/home/f002s7f/Desktop/drives/54054-Linux/54054dua/idata/cbobak/"
outdir_2<-"/dsrives/drive1/home/f002s7f/Desktop/drives/54054-Linux/54054dua/idata/cbobak/participation_analysis/"
outdir_spillover<-paste(outdir,"spillover_analysis/",sep="")

setwd(outdir_spillover)
stepwise_network_list<-readRDS("phys_network_list.RData")
npi_trial_data<-read.csv(paste(outdir,"npi_trial_data.csv",sep=""))
npi_trial_data<-npi_trial_data %>% rename(HospitalID = "trial_hospitalid")

# 1) Ensure consistent vertex order across steps (union of vertices)
all_npi <- Reduce(union, lapply(stepwise_network_list, function(g) V(g)$name))

make_W <- function(g, all_ids) {
  # adjacency as sparse (unweighted OR weighted—your choice; we’ll row-normalize)
  A <- as_adj(g, attr = "weight", sparse = TRUE)  # if you want weights
  # align to full set with zero rows/cols for missing NPIs
  cur <- rownames(A)
  idx <- match(all_ids, cur)
  # build aligned matrix
  A_full <- Matrix(0, length(all_ids), length(all_ids), sparse = TRUE)
  rownames(A_full) <- colnames(A_full) <- all_ids
  if (!all(is.na(idx))) {
    A_full[cbind(match(cur, all_ids), match(cur, all_ids))] <- 0  # no-op; keeps diag 0
    A_full[match(cur, all_ids), match(cur, all_ids)] <- A
  }
  # zero diag (safety)
  diag(A_full) <- 0
  # row-stochastic with isolated handling
  rs <- Matrix::rowSums(A_full)
  iso <- which(rs == 0)
  if (length(iso)) A_full[cbind(iso, iso)] <- 1
  A_full <- A_full / Matrix::rowSums(A_full)
  A_full
}

W_list <- lapply(stepwise_network_list, make_W, all_ids = all_npi)

# 2) Metadata vectors (named by NPI)
npi <- all_npi
trial_step <- setNames(npi_trial_data$trial_step[match(npi, npi_trial_data$npi)], npi)
trial_step[is.na(trial_step)]<-Inf #avoids errors later
hospital   <- setNames(all_npi_hosp_match$HospitalID[match(npi, all_npi_hosp_match$NPI)], npi)

# 3) Bundle + save
network_panel <- list(W_list = W_list, npi = npi, trial_step = trial_step, hospital = hospital)
saveRDS(network_panel, file = "network_panel.rds")
