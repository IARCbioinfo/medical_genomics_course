library(MOFA2)
file <- system.file("extdata", "test_data.RData", package = "MOFA2")
load(file)
# Create the MOFA object
MOFAmodel <- create_mofa(dt)
# Prepare the MOFA object with default options
MOFAmodel <- prepare_mofa(MOFAmodel)
# Run the MOFA model
MOFAmodel <- run_mofa(MOFAmodel)
