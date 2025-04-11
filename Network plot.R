
# Load necessary libraries
library(tidyverse)
library(igraph)

# Load your datasets
d <- file.choose()
d <- read.csv(d, header = TRUE)
p <- file.choose()   
p <- read.csv(p, header = TRUE)

# Load necessary libraries
library(igraph)

# Prepare the proximity data by removing the first column (which is "X")
proximity_data <- p[, -1]

# Create an empty data frame for storing the edge list
edges <- data.frame(from = character(0), to = character(0), weight = numeric(0))

# Loop through the upper triangle of the matrix to avoid duplicate edges
for (i in 1:(ncol(proximity_data) - 1)) {
  for (j in (i + 1):ncol(proximity_data)) {
    if (proximity_data[i, j] > 0) {  # Only consider non-zero proximity
      edges <- rbind(edges, data.frame(from = colnames(proximity_data)[i], 
                                       to = colnames(proximity_data)[j], 
                                       weight = proximity_data[i, j]))
    }
  }
}

# Create graph from the edge list
g <- graph_from_data_frame(edges, directed = FALSE)

# Check the first few rows of the dataset and symptomatic data
print(head(d)) 

symptomatic <- d$symptomatic  # Adjust the column name if necessary

# Ensure vertex names in the graph match the identities from the dataset
vertex_names <- V(g)$name

# Check the first few vertex names
print(head(vertex_names))

# If necessary, reorder 'symptomatic' to match the vertex order
# We ensure that the 'symptomatic' values match the correct individuals
symptomatic_ordered <- symptomatic[match(vertex_names, d$identity)]

# Assign colors based on the symptomatic status
V(g)$color <- ifelse(symptomatic_ordered == 1, "red", "green")

# Ensure that the correct layout (positioning) is used
layout_matrix <- layout_with_fr(g)

# Visualize the network with correct node positions and colors
plot(g, 
     layout = layout_matrix,  # Apply Fruchterman-Reingold layout for positioning
     vertex.size = 15, 
     vertex.label.cex = 0.7, 
     edge.width = E(g)$weight * 10,  # Scale edge width by proximity
     vertex.color = V(g)$color,     # Apply the color from the 'symptomatic' data
     main = "Proximity Network with Symptomatic Status")



#Maybe relevant for counting observations if dealing with raw data


# Load required packages
library(dplyr)

# Step 1: Count the number of observations per chimpanzee (based on 'identity')
d_observations <- d |>
  group_by(identity) |>
  summarise(observations = n())

# Merge the observations back into the original dataset 'd'
d <- left_join(d, d_observations, by = "identity")

# Step 2: Linear Regression for Centrality vs. Number of Observations
lm_centrality <- lm(strength.centrality ~ observations, data = d)
summary(lm_centrality)

# Step 3: Linear Regression for Symptomatic vs. Number of Observations
lm_symptomatic <- lm(symptomatic ~ observations, data = d)
summary(lm_symptomatic)



var(d$observations)
cor(d$observations, d$strength.centrality)
# Try removing the `observations` variable and fit the model
lm_symptomatic_simple <- lm(symptomatic ~ 1, data = d)  # Only intercept
summary(lm_symptomatic_simple)






# Step 1: Prepare the proximity data by removing the first column (which is "X")
proximity_data <- p[, -1]

# Create an empty data frame for storing the edge list
edges <- data.frame(from = character(0), to = character(0), weight = numeric(0))

# Loop through the upper triangle of the matrix to avoid duplicate edges
for (i in 1:(ncol(proximity_data) - 1)) {
  for (j in (i + 1):ncol(proximity_data)) {
    if (proximity_data[i, j] > 0) {  # Only consider non-zero proximity
      edges <- rbind(edges, data.frame(from = colnames(proximity_data)[i], 
                                       to = colnames(proximity_data)[j], 
                                       weight = proximity_data[i, j]))
    }
  }
}

# Create graph from the edge list
g <- graph_from_data_frame(edges, directed = FALSE)

# Ensure vertex names in the graph match the identities from the dataset
vertex_names <- V(g)$name

# Step 2: Check if all 'from' and 'to' nodes exist in the graph
node_pairs <- combn(V(g)$name, 2, simplify = FALSE)

# Debugging: Check if the nodes in the node pairs exist in the graph
invalid_nodes <- unlist(node_pairs)[!unlist(node_pairs) %in% vertex_names]
if (length(invalid_nodes) > 0) {
  cat("Invalid nodes found in node pairs: ", unique(invalid_nodes), "\n")
}

# Ensure no NA values are passed to the function
node_pairs <- node_pairs[!sapply(node_pairs, function(pair) any(is.na(pair)))]

# Function to calculate k-path length between two nodes
k_path_length <- function(graph, node1, node2, k) {
  # Ensure that node1 and node2 are valid vertex names in the graph
  if (!(node1 %in% V(graph)$name) | !(node2 %in% V(graph)$name)) {
    cat("Invalid vertex names:", node1, node2, "\n")
    return(NA)
  }
  
  paths <- all_shortest_paths(graph, from = node1, to = node2, weights = E(graph)$weight)
  
  # Filter paths by length
  valid_paths <- Filter(function(x) length(x) == k + 1, paths$res)
  return(length(valid_paths))  # Return the number of k-paths
}

# Step 3: Calculate k-paths for the observed network
k <- 3  # Set the k-path length

# Debugging: Check the first few node pairs
print(head(node_pairs))

# Calculate observed k-path lengths
observed_k_paths <- sapply(node_pairs, function(pair) {
  k_path_length(g, pair[1], pair[2], k)
})

# Step 4: Generate random networks for null hypothesis testing
# Function to generate random network by preserving degree distribution
generate_random_network <- function(g) {
  sample_degseq(degree(g), method = "vl")
}

# Step 5: Compute k-path lengths for random networks
n_random <- 1000  # Number of random networks to generate
random_k_paths <- replicate(n_random, {
  random_g <- generate_random_network(g)
  random_k_paths <- sapply(node_pairs, function(pair) {
    k_path_length(random_g, pair[1], pair[2], k)
  })
  return(mean(random_k_paths))  # Mean k-path length for the random network
})

# Step 6: Compare the observed k-path length to the distribution of random k-path lengths
observed_mean_k_path <- mean(observed_k_paths)  # Mean k-path length in the observed network
p_value <- mean(random_k_paths >= observed_mean_k_path)  # p-value for the significance test

# Step 7: Print the result
cat("Observed mean k-path length:", observed_mean_k_path, "\n")
cat("p-value for k-path test:", p_value, "\n")



summary(random_graph)  # Check for basic properties of the random graph
summary(random_k_paths)  # Check the distribution of random k-path lengths
# p-value is the fraction of random graph k-paths that are greater than or equal to the observed value
p_value <- mean(random_k_paths >= observed_value)
cat("p-value for k-path test:", p_value, "\n")
# Remove NA values from random_k_paths before calculating the p-value
random_k_paths_clean <- na.omit(random_k_paths)
p_value <- mean(random_k_paths_clean >= observed_value)
cat("p-value for k-path test:", p_value, "\n")



print(V(g)$name)  # Check available vertex names
node_pairs <- combn(na.omit(V(g)$name), 2, simplify = FALSE)
vertex_names <- V(g)$name
valid_nodes <- na.omit(match(vertex_names, d$identity))


print(setdiff(vertex_names, d$identity))  # Names in graph but NOT in dataset
print(setdiff(d$identity, vertex_names))  # Names in dataset but NOT in graph
d$identity <- as.character(d$identity)  # Ensure character format
valid_nodes <- match(vertex_names, d$identity)
valid_nodes <- na.omit(valid_nodes)  # Remove NAs
d$identity <- trimws(d$identity)  # Remove leading/trailing spaces
vertex_names <- trimws(vertex_names)


print(sum(is.na(valid_nodes)))  # Count NA values after match
print(sum(is.na(d$identity)))  # Check if d$identity had NAs
print(sum(is.na(vertex_names)))  # Check if vertex_names had NAs
valid_nodes <- match(vertex_names, d$identity, nomatch = 0)  # Returns 0 for unmatched values
valid_nodes <- valid_nodes[valid_nodes > 0]  # Remove zeros (previously unmatched)
print(valid_nodes)  # Check the matched indices



valid_vertex_set <- V(g)[valid_nodes]
print(valid_vertex_set)  # Confirm valid node selection
valid_vertex_names <- vertex_names[valid_nodes]
print(valid_vertex_names)  # Check the extracted names
valid_vertex_set <- V(g)[name %in% valid_vertex_names]
print(valid_vertex_set)
print(from)
print(class(from))



as_igraph_vs(graph, from)

valid_vertex_set <- V(g)[name %in% valid_vertex_names]
print(valid_vertex_set)




# Generate a random graph with the same degree sequence using sample_degseq
deg_seq <- degree(g)  # Get the degree sequence from the actual graph
random_graph_degseq <- sample_degseq(deg_seq)  # Sample a random graph with the same degree sequence

# Function to compute the mean weighted path (adjust if needed)
mean_weighted_path <- function(graph, valid_nodes) {
  edge_weights <- E(graph)$weight
  return(mean(edge_weights))  # Modify this function if you need specific k-paths
}

# Calculate the observed k-path length for the actual graph
observed_value <- mean_weighted_path(g, valid_vertex_set)

# Randomize the graph 1000 times and calculate k-path lengths for each randomization
random_k_paths <- replicate(1000, {
  # Re-generate the random graph with the same degree sequence
  random_graph_degseq <- sample_degseq(deg_seq)
  
  # Calculate the k-path length for the random graph
  random_k_path <- mean_weighted_path(random_graph_degseq, valid_vertex_set)
  return(random_k_path)
})

# Clean up random k-paths by removing NAs if any are present
random_k_paths_clean <- na.omit(random_k_paths)

# Calculate p-value: Fraction of random k-path lengths greater than or equal to the observed value
p_value <- mean(random_k_paths_clean >= observed_value)

# Output the observed k-path length and p-value
cat("Observed mean k-path length:", observed_value, "\n")
cat("p-value for k-path test:", p_value, "\n")

# Handle edge weights and missing values (if needed)
E(g)$weight <- as.numeric(E(g)$weight)
sum(is.na(E(g)$weight))  # Check for any NAs
E(g)$weight[is.na(E(g)$weight)] <- 0  # Replace NAs with 0, or use another imputation method





# Check for warnings
warnings()

# Print the cleaned random k-paths
print(random_k_paths_clean)

# Modified p-value calculation with a check for empty random k-paths
if (length(random_k_paths_clean) > 0) {
  p_value <- mean(random_k_paths_clean >= observed_value)
} else {
  p_value <- NA  # Assign a default value if no valid random k-paths
}

# Output the observed k-path length and p-value
cat("Observed mean k-path length:", observed_value, "\n")
cat("p-value for k-path test:", p_value, "\n")

# Check for edge weights and missing values
E(g)$weight <- as.numeric(E(g)$weight)
sum(is.na(E(g)$weight))  # Check for any missing edge weights
E(g)$weight[is.na(E(g)$weight)] <- 0  # Replace NAs with 0, or use another imputation method




# Load necessary libraries
library(igraph)

# Assuming 'random_k_paths' is your data containing the paths
# Step 1: Inspect the structure of 'random_k_paths' to understand its type
print(str(random_k_paths))

# Step 2: Try to convert 'random_k_paths' to numeric properly
# If it's a character or factor vector, try converting it explicitly
random_k_paths <- as.numeric(as.character(random_k_paths))

# Step 3: Check if the conversion worked and if there are any NAs
print(sum(is.na(random_k_paths)))  # Check how many NA values
print(random_k_paths)  # Inspect the first few elements

# Step 4: Remove NA values from 'random_k_paths'
random_k_paths_clean <- na.omit(random_k_paths)

# Step 5: Ensure that the cleaned data is numeric and not empty
if (length(random_k_paths_clean) == 0) {
  stop("No valid numeric data available after cleaning.")
} else {
  print("Data cleaned successfully!")
}

# Step 6: Calculate the mean of the cleaned random_k_paths
mean_random_k_path <- mean(random_k_paths_clean)
print(paste("Mean of random k-paths: ", mean_random_k_path))

# Optional: If working with igraph object and needing edge weights
# If 'graph' is your igraph object, get the edge weights as follows
# Replace 'graph' with your actual graph object
edge_weights <- E(graph)$weight

# Check if edge weights are valid
edge_weights <- as.numeric(edge_weights)

if (length(edge_weights) > 0 && !all(is.na(edge_weights))) {
  mean_edge_weight <- mean(edge_weights, na.rm = TRUE)
  print(paste("Mean of edge weights: ", mean_edge_weight))
} else {
  print("No valid edge weights found.")
}




#Descriptive stats

# Descriptive statistics for centrality metrics in the dataset 'd'
summary_stats <- data.frame(
  Metric = c("Strength", "Eigenvector Centrality", "Flow Betweenness"),
  Mean = c(mean(d$strength.centrality, na.rm=TRUE), 
           mean(d$eigenvector.centrality, na.rm=TRUE), 
           mean(d$flow.betweenness, na.rm=TRUE)),
  Median = c(median(d$strength.centrality, na.rm=TRUE), 
             median(d$eigenvector.centrality, na.rm=TRUE), 
             median(d$flow.betweenness, na.rm=TRUE)),
  SD = c(sd(d$strength.centrality, na.rm=TRUE), 
         sd(d$eigenvector.centrality, na.rm=TRUE), 
         sd(d$flow.betweenness, na.rm=TRUE)),
  Min = c(min(d$strength.centrality, na.rm=TRUE), 
          min(d$eigenvector.centrality, na.rm=TRUE), 
          min(d$flow.betweenness, na.rm=TRUE)),
  Max = c(max(d$strength.centrality, na.rm=TRUE), 
          max(d$eigenvector.centrality, na.rm=TRUE), 
          max(d$flow.betweenness, na.rm=TRUE)),
  Range = c(range(d$strength.centrality, na.rm=TRUE)[2] - range(d$strength.centrality, na.rm=TRUE)[1],
            range(d$eigenvector.centrality, na.rm=TRUE)[2] - range(d$eigenvector.centrality, na.rm=TRUE)[1],
            range(d$flow.betweenness, na.rm=TRUE)[2] - range(d$flow.betweenness, na.rm=TRUE)[1])
)

# Print the descriptive statistics
print(summary_stats)

# Histogram of Strength Centrality
hist(d$strength.centrality, 
     xlab="Strength (Degree Centrality)", 
     main="Histogram of Strength", 
     col="lightblue", 
     border="black")




#Boxplots

# Load the necessary library
library(ggplot2)

# Prepare the data for the boxplot
centrality_data <- data.frame(
  Strength = d$strength.centrality,
  Eigenvector = d$eigenvector.centrality,
  FlowBetweenness = d$flow.betweenness
)

# Create the boxplot for Strength, Eigenvector Centrality, and Flow Betweenness
ggplot(centrality_data, aes(x = factor(1), y = Strength)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(x = "Strength (Degree Centrality)", y = "Centrality Value") +
  theme_minimal() +
  ggtitle("Boxplot of Strength Centrality")

# Boxplot for Eigenvector Centrality
ggplot(centrality_data, aes(x = factor(1), y = Eigenvector)) +
  geom_boxplot(fill = "lightgreen", color = "black") +
  labs(x = "Eigenvector Centrality", y = "Centrality Value") +
  theme_minimal() +
  ggtitle("Boxplot of Eigenvector Centrality")

# Boxplot for Flow Betweenness
ggplot(centrality_data, aes(x = factor(1), y = FlowBetweenness)) +
  geom_boxplot(fill = "lightcoral", color = "black") +
  labs(x = "Flow Betweenness", y = "Centrality Value") +
  theme_minimal() +
  ggtitle("Boxplot of Flow Betweenness")



head(d)

# Load necessary package
library(car)  # for VIF

# Fit logistic regression
model <- glm(symptomatic ~ strength.centrality + dominance.score + age, 
             family = binomial, data = d)

# Display model summary
summary(model)

# Load necessary package
library(boot)

# Function to shuffle symptomatic variable and refit model
permute_model <- function(data, indices) {
  permuted_data <- data
  permuted_data$symptomatic <- sample(data$symptomatic)  # Shuffle labels
  perm_model <- glm(symptomatic ~ strength.centrality + dominance.score + age, 
                    family = binomial, data = permuted_data)
  return(coef(perm_model))  # Store coefficients
}

# Run permutation test with 30,000 iterations
set.seed(123)
perm_results <- boot(data = d, statistic = permute_model, R = 30000)

# Extract permutation-based p-values
perm_p_values <- apply(perm_results$t, 2, function(x) mean(abs(x) >= abs(coef(model))))

# Display permutation-based p-values
perm_p_values

#running vif
library(car)
vif(glm(symptomatic ~ strength.centrality + dominance.score + age, 
        family = binomial, data = d))


# Permutation-based regression for Dominance Rank on Strength Centrality
set.seed(123)  # for reproducibility
perm_dominance <- replicate(30000, {
  d$permuted_outcome <- sample(d$strength.centrality)
  fit_perm <- glm(permuted_outcome ~ dominance.score, family = binomial, data = d)
  coef(fit_perm)[2]  # Get the coefficient for dominance.score
})

# Calculate p-value for Dominance Score
obs_dominance <- coef(glm(strength.centrality ~ dominance.score, family = binomial, data = d))[2]
perm_p_value_dominance <- mean(abs(perm_dominance) >= abs(obs_dominance))

# Permutation-based regression for Age on Strength Centrality
perm_age <- replicate(30000, {
  d$permuted_outcome <- sample(d$strength.centrality)
  fit_perm <- glm(permuted_outcome ~ age, family = binomial, data = d)
  coef(fit_perm)[2]  # Get the coefficient for age
})

# Calculate p-value for Age
obs_age <- coef(glm(strength.centrality ~ age, family = binomial, data = d))[2]
perm_p_value_age <- mean(abs(perm_age) >= abs(obs_age))

# Print results
print(paste("Permutation p-value for Dominance Rank: ", perm_p_value_dominance))
print(paste("Permutation p-value for Age: ", perm_p_value_age))


#Answer: Permutation p-value for Dominance Rank: 0This indicates that the effect of dominance rank on strength centrality is highly significant and unlikely to have occurred by chance. The very low p-value suggests a strong relationship between dominance rank and strength centrality in the social network.Permutation p-value for Age: 4e-04Similarly, the effect of age on strength centrality is also significant. The p-value is quite small, suggesting that age does play a role in determining strength centrality within the network, albeit perhaps not as strongly as dominance rank.

#Interpretation: These results imply that both dominance rank and age are important predictors of strength centrality in the network. The dominance rank, in particular, has a very strong association, and age also significantly contributes to the centrality score. These findings suggest that higher dominance rank and certain age-related factors (perhaps social experience or physiological changes) might affect a chimpanzee's position and influence within the network.





colnames(d)

# Create ordinal rank based on dominance score
d$ordinal.rank <- rank(d$dominance.score, ties.method = "average")

# Logistic regression with glm function
model <- glm(symptomatic ~ strength.centrality + ordinal.rank + age, 
             family = binomial, data = d)

# Display model summary
summary(model)

# Assessing collinearity using VIF
library(car)
vif(model)

# Permutation p-values for variable importance
# Function for permutation
permute_glm <- function(model, data, n_iter = 30000) {
  # Store original regression coefficients
  original_coefs <- coef(model)
  
  # Set up a vector to store permutation coefficients
  permuted_coefs <- matrix(NA, nrow = n_iter, ncol = length(original_coefs))
  
  # Permute the response variable and refit the model
  for (i in 1:n_iter) {
    permuted_data <- data
    permuted_data$symptomatic <- sample(permuted_data$symptomatic)
    
    permuted_model <- glm(symptomatic ~ strength.centrality + ordinal.rank + age, 
                          family = binomial, data = permuted_data)
    
    permuted_coefs[i, ] <- coef(permuted_model)
  }
  
  # Calculate p-values based on permutation results
  p_values <- apply(permuted_coefs, 2, function(x) {
    mean(abs(x) >= abs(original_coefs))
  })
  
  return(p_values)
}

# Get permutation p-values
perm_p_values <- permute_glm(model, d)
print(perm_p_values)





# Visualization: Boxplot for Strength Centrality
ggplot(d, aes(x = factor(1), y = strength.centrality)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(x = "Strength Centrality", y = "Centrality Value") +
  theme_minimal() +
  ggtitle("Boxplot of Strength Centrality")

# Visualization: Boxplot for Eigenvector Centrality
ggplot(d, aes(x = factor(1), y = eigenvector.centrality)) +
  geom_boxplot(fill = "lightgreen", color = "black") +
  labs(x = "Eigenvector Centrality", y = "Centrality Value") +
  theme_minimal() +
  ggtitle("Boxplot of Eigenvector Centrality")

# Visualization: Boxplot for Flow Betweenness
ggplot(d, aes(x = factor(1), y = flow.betweenness)) +
  geom_boxplot(fill = "lightcoral", color = "black") +
  labs(x = "Flow Betweenness", y = "Centrality Value") +
  theme_minimal() +
  ggtitle("Boxplot of Flow Betweenness")








#Raw Data
edge_list <- p_raw %>%
  select(focalID, prox5) %>%
  rename(partnerID = prox5) 


edge_list <- as.data.frame(edge_list)


head(edge_list)

network <- graph_from_data_frame(edge_list, directed = FALSE)

degree_centrality <- degree(network)


print(degree_centrality)


# Scale the node size
node_size <- degree_centrality / max(degree_centrality) * 10  

# Visualize network
plot(network, vertex.size=node_size, 
     vertex.label=V(network)$name, 
     vertex.label.cex=0.8, 
     edge.width=1.5, 
     main="Social Network Based on Proximity",
     layout=layout_with_fr)  


flow_betweenness_vals <- betweenness(network, directed = TRUE, normalized = TRUE)
print("Flow Betweenness Centrality Values:")
print(flow_betweenness_vals)

eigenvector_vals <- eigen_centrality(network, directed = TRUE)$vector
print("Eigenvector Centrality Values:")
print(eigenvector_vals)







  


