# install_deps.R
# Check and install required R packages for classprediction pipeline

pkgs <- c(
  "dplyr",
  "stringr",
  "reshape2",
  "ggplot2",
  "igraph",
  "ggraph",
  "gridExtra",
  "scales",
  "data.table",
  "tidyr",
  "rpart",
  "purrr"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing missing package: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  } else {
    message("Package already installed: ", p)
  }
}

message("\nAll dependencies are installed and available.\n")

# Print versions for reproducibility
for (p in pkgs) {
  cat(p, "version:", as.character(packageVersion(p)), "\n")
}
