num_samples = 330_000  # Total number of samples
num_features = 80_000  # Example value for 'n', the number of float values per sample
float_size_bytes = 8  # Size of a float in bytes (64-bit float)

# Calculate the total memory usage in bytes
total_memory_bytes = num_samples * num_features * float_size_bytes

# Convert bytes to gigabytes for easier interpretation
total_memory_gb = total_memory_bytes / (1024 ** 3)

print(f'Total memory usage: {total_memory_gb:.2f} GB')

# Estimate linear regression memory usage

# Calculate the size of the X^T X matrix
xtx_matrix_size = num_features * num_features

# Calculate the memory usage of the X^T X matrix in bytes
xtx_memory_bytes = xtx_matrix_size * float_size_bytes

# Convert bytes to gigabytes for easier interpretation
xtx_memory_gb = xtx_memory_bytes / (1024 ** 3)

print(f'Memory usage for X^T X matrix: {xtx_memory_gb:.2f} GB')