# Load necessary libraries
library(dplyr)    # For data manipulation
library(Matrix)   # For handling sparse matrices
library(igraph)   # For network construction and analysis

# Build data frame containing NPIs involved in the ACP intervention
acp_df <- data.frame(acp_npi)

# Filter encounter data for physicians participating in the trial within the trial period
npi_trial_location <- encounters %>%
  filter(between(DateOfService, as.Date("2020-10-06"), as.Date("2021-03-01")),
         NPI %in% acp_npi) %>%
  group_by(HospitalID, NPI) %>%
  summarize(numEncounters = n()) %>%
  group_by(NPI) %>%
  filter(numEncounters == max(numEncounters)) %>%
  ungroup() %>%
  select(NPI, HospitalID)

# Ensure DateOfService is in date format
encounters$DateOfService <- as.Date(encounters$DateOfService)

# Filter encounter data to include relevant NPIs and timeframe for the study
encounters_network <- encounters %>%
  filter(NPI %in% acp_npi,
         between(DateOfService, as.Date("2019-01-01"), as.Date("2021-06-01")))

# Create the initial sparse adjacency matrix of NPIs by patient stay
m <- table(encounters_network$NPI, encounters_network$patient_stay_id)

# Optionally divide into manageable chunks if the matrix is too large for memory
# Here we divide the matrix into chunks of 100,000 patient stays
m1 <- m[, 1:100000]
m2 <- m[, 100001:200000]
m3 <- m[, 200001:300000]
m4 <- m[, 300001:400000]
m5 <- m[, 400001:500000]
m6 <- m[, 500001:600000]
m7 <- m[, 600001:700000]
m8 <- m[, 700001:ncol(m)]

# List of all matrix chunks for iterative processing
mlist <- list(m1, m2, m3, m4, m5, m6, m7, m8)

# Initialize an overlap matrix to store cumulative shared patient counts
for (i in 1:length(mlist)) {
  m <- as(mlist[[i]], "sparseMatrix")  # Convert to sparse matrix format
  
  # On the first iteration, initialize the overlap matrix
  if (i == 1) {
    overlapMatrix <- m %*% t(m)
  } else {
    overlapMatrix <- overlapMatrix + (m %*% t(m))  # Accumulate overlap
  }
}

# Set row and column names for the overlap matrix to NPIs
colnames(overlapMatrix) <- rownames(overlapMatrix) <- rownames(m)

# Create a graph object from the overlap matrix (undirected and weighted)
g <- graph_from_adjacency_matrix(overlapMatrix, mode = "undirected", diag = FALSE, weighted = TRUE)

# Generate coordinates for the network layout using Graphopt layout algorithm
coords <- layout_with_graphopt(g, niter = 500, charge = 0.01, mass = 30, spring.length = 0, spring.constant = 1)

# Assign row names of the coordinates to match the graph vertices
row.names(coords) <- V(g)$name

# Plot the graph (this part can be omitted or refined for large datasets)
plot(g, vertex.size = 7, vertex.label.cex = 0.05, layout = coords)

# Define the start and end dates of network snapshots for each trial step
StartNetworks <- as.Date("2020-04-06")
EndNetworks <- as.Date(c("2020-10-05", "2020-10-30", "2020-12-05", "2020-12-30", "2021-01-26", "2021-03-01"))

# Function to divide the adjacency matrix into smaller chunks for processing
divide_m <- function(m, chunk_size = 100000) {
  nc <- ncol(m)
  num_chunks <- ceiling(nc / chunk_size)
  
  mlist <- vector("list", num_chunks)
  
  for (i in 1:num_chunks) {
    start <- (i - 1) * chunk_size + 1
    end <- min(i * chunk_size, nc)
    mlist[[i]] <- m[, start:end]
  }
  
  return(mlist)
}

# Function to create a physician interaction graph for a given time period
makePhysGraph <- function(df, start, end) {
  
  # Filter encounters for the specified time period
  df <- df %>% filter(between(DateOfService, start, end))
  
  # Create the NPI by patient stay matrix
  m <- table(df$NPI, df$patient_stay_id)
  
  # Divide the matrix into chunks for processing
  mlist <- divide_m(m)
  
  # Initialize the overlap matrix
  for (i in 1:length(mlist)) {
    m <- as(mlist[[i]], "sparseMatrix")
    
    if (i == 1) {
      overlapMatrix <- m %*% t(m)
    } else {
      overlapMatrix <- overlapMatrix + (m %*% t(m))
    }
  }
  
  # Create the interaction graph from the overlap matrix
  g <- graph_from_adjacency_matrix(overlapMatrix, mode = "undirected", diag = FALSE, weighted = TRUE)
  
  return(g)
}

# Initialize lists to store network snapshots and graph statistics
stepwise_network_list <- list()
step <- 0

# Vectors to store network statistics across time points
densityV <- c()
transitV <- c()
diameterV <- c()
min_deg <- c()
med_deg <- c()
max_deg <- c()
sd_deg <- c()

# Loop through each step in the trial and compute network statistics
for (i in 1:length(EndNetworks)) {
  # Create a network for the current time window
  g <- makePhysGraph(encounters_network, StartNetworks, EndNetworks[i])
  
  # Store the network in the list of stepwise networks
  stepwise_network_list[[paste("Step", step)]] <- g
  
  # Calculate and store network statistics
  densityV[step + 1] <- edge_density(g, loops = FALSE)
  transitV[step + 1] <- transitivity(g, type = "global")
  diameterV[step + 1] <- diameter(g, directed = FALSE)
  deg <- degree(g, mode = "all")
  min_deg[step + 1] <- min(deg)
  max_deg[step + 1] <- max(deg)
  med_deg[step + 1] <- median(deg)
  sd_deg[step + 1] <- sd(deg)
  
  step <- step + 1
}

# Summarize network statistics in a data frame
summary_stats <- data.frame(
  Step = 0:5,
  Density = densityV,
  Transitivity = transitV,
  Diameter = diameterV,
  Minimum_Degree = min_deg,
  Median_Degree = med_deg,
  Maximum_Degree = max_deg,
  Degree_Std_Dev = sd_deg
)

# Create a data frame to track physician trial effects across time steps
phys_trial_effects <- data.frame(NPI = V(stepwise_network_list[[6]])$name)
rownames(npi_trial_data) <- npi_trial_data$npi

# Calculate cumulative exposure to the intervention for each physician across steps
for (step in 1:5) {
  W <- as_adj(stepwise_network_list[[step + 1]])  # Get adjacency matrix for the step
  Post <- ifelse(npi_trial_data[W@Dimnames[[1]], "trial_step"] <= step, 1, 0)
  Post[is.na(Post)] <- 0
  WxPost <- (W %*% Post)[, 1]
  WxPost <- WxPost[phys_trial_effects$NPI]
  
  # Add the calculated WxPost values to the data frame
  phys_trial_effects <- within(phys_trial_effects, assign(paste("WxPost", step, sep = ""), WxPost))
}

# Calculate row-wise stochastic exposure (strict prior exposure modeling)
phys_trial_effects_rw <- data.frame(NPI = V(stepwise_network_list[[6]])$name)

for (step in 1:5) {
  W <- as_adj(stepwise_network_list[[step]])
  rs <- rowSums(W)
  
  # Ensure no rows have zero sum (i.e., disconnected nodes)
  if (any(rs == 0)) {
    rs[rs==0]<-1
  }
  
  # Normalize rows to ensure stochastic interpretation
  W <- W /rs
  
  # Calculate exposure using trial step data
  Post <- ifelse(npi_trial_data[W@Dimnames[[1]], "trial_step"] < step, 1, 0)
  Post[is.na(Post)] <- 0
  WxPost <- (W %*% Post)[, 1]
  WxPost <- WxPost[phys_trial_effects_rw$NPI]
  
  # Add the calculated stochastic WxPost values to the data frame
  phys_trial_effects_rw <- within(phys_trial_effects_rw, assign(paste("WxPost", step, sep = ""), WxPost))
}
