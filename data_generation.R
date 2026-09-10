library(bnlearn)

# -------------------------------------------------
# ALARM NETWORK (GROUND-TRUTH)
# Sub-critical sparse benchmark (p = 37): seed 123
# -------------------------------------------------

# Load the downloaded alarm data
alarm_true_bn <- readRDS("alarm.rds")
print(alarm_true_bn)

# Distinct seeds per network to avoid any shared-stream correlation across benchmarks
set.seed(123)

# Generate 5 datasets
data_500 <- rbn(alarm_true_bn, n=500)
data_1000 <- rbn(alarm_true_bn, n=1000)
data_2500 <- rbn(alarm_true_bn, n=2500)
data_5000 <- rbn(alarm_true_bn, n=5000)
data_10000 <- rbn(alarm_true_bn, n=10000)

# Save as CSV files
write.csv(data_500, "data/alarm_500.csv", row.names = FALSE)
write.csv(data_1000, "data/alarm_1000.csv", row.names = FALSE)
write.csv(data_2500, "data/alarm_2500.csv", row.names = FALSE)
write.csv(data_5000, "data/alarm_5000.csv", row.names = FALSE)
write.csv(data_10000, "data/alarm_10000.csv", row.names = FALSE)


# ---------------------------------------------------
# HEPAR II NETWORK
# Super-critical dense benchmark (p = 70): seed 234
# ---------------------------------------------------

# Load the downloaded hepar2 data
hepar_true_bn <- readRDS("hepar2.rds")
print(hepar_true_bn)

set.seed(234)

# Generate 5 datasets
hepar_data_500 <- rbn(hepar_true_bn, n=500)
hepar_data_1000 <- rbn(hepar_true_bn, n=1000)
hepar_data_2500 <- rbn(hepar_true_bn, n=2500)
hepar_data_5000 <- rbn(hepar_true_bn, n=5000)
hepar_data_10000 <- rbn(hepar_true_bn, n=10000)

# Save as CSV files
write.csv(hepar_data_500, "data/hepar_500.csv", row.names = FALSE)
write.csv(hepar_data_1000, "data/hepar_1000.csv", row.names = FALSE)
write.csv(hepar_data_2500, "data/hepar_2500.csv", row.names = FALSE)
write.csv(hepar_data_5000, "data/hepar_5000.csv", row.names = FALSE)
write.csv(hepar_data_10000, "data/hepar_10000.csv", row.names = FALSE) 

# --------------------------------------------------
# HAILFINDER NETWORK
# Intermediate density benchmark (p = 56): seed 345
# --------------------------------------------------

# Load the downloaded hailfinder data
hailfinder_true_bn <- readRDS("hailfinder.rds")
print(hailfinder_true_bn)

set.seed(345)

# Generate 5 datasets
hailfinder_data_500 <- rbn(hailfinder_true_bn, n=500)
hailfinder_data_1000 <- rbn(hailfinder_true_bn, n=1000)
hailfinder_data_2500 <- rbn(hailfinder_true_bn, n=2500)
hailfinder_data_5000 <- rbn(hailfinder_true_bn, n=5000)
hailfinder_data_10000 <- rbn(hailfinder_true_bn, n=10000)

# Save as CSV files
write.csv(hailfinder_data_500, "data/hailfinder_500.csv", row.names = FALSE)
write.csv(hailfinder_data_1000, "data/hailfinder_1000.csv", row.names = FALSE)
write.csv(hailfinder_data_2500, "data/hailfinder_2500.csv", row.names = FALSE)
write.csv(hailfinder_data_5000, "data/hailfinder_5000.csv", row.names = FALSE)
write.csv(hailfinder_data_10000, "data/hailfinder_10000.csv", row.names = FALSE) 

print("Week 1 Data Generation Complete! All CSV files are saved safely.")

