import cv2
import numpy as np
import matplotlib.pyplot as plt
import os

# Function to calculate the mean pixel value in the specified ROI
def calculate_mean_pixel_value_in_roi(image, roi_coords):
    x1, y1, x2, y2 = roi_coords
    # Extract the ROI from the image
    roi = image[y1:y2, x1:x2]
    # Convert the ROI to grayscale
    gray_roi = cv2.cvtColor(roi, cv2.COLOR_BGR2GRAY)
    # Calculate the mean of the grayscale ROI
    mean_value = np.mean(gray_roi)
    return mean_value

# Directory where your images are stored
image_directory ="path/to/spectralon/images"

# List your images in the directory (ensure they are sorted correctly by exposure time)
image_files = sorted([os.path.join(image_directory, f) for f in os.listdir(image_directory) if f.endswith('.tif')])

# Define the coordinates of the ROI (example: top-left corner (x1, y1) and bottom-right corner (x2, y2))
roi_coords = (600, 480, 770, 570)  # Replace with the actual coordinates for your white standard

# Placeholder for storing the mean pixel values and exposure times
mean_values = []
exposure_times = []

# Iterate over the images and process each one
for image_file in image_files:
    # Read the image
    image = cv2.imread(image_file)

    # Calculate the mean pixel value in the ROI
    mean_value = calculate_mean_pixel_value_in_roi(image, roi_coords)
    mean_values.append(mean_value)

    # Extract exposure time from the filename (assuming filename format contains exposure time information)
    # You may need to adjust this part depending on how exposure times are encoded in your filenames
    exposure_time = int(os.path.splitext(os.path.basename(image_file))[0].split('_')[-1][:-3])  # Example extraction
    exposure_times.append(exposure_time)

# Plotting the results
plt.figure(figsize=(10, 6))
plt.scatter(exposure_times, mean_values, marker='o', linestyle='-', color='b')
plt.xlabel('Exposure Time (ms)')
plt.ylabel('Mean Pixel Value in ROI')
plt.title('Mean Pixel Value in ROI vs. Exposure Time')
plt.grid(True)
plt.show()


from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

# Sample 1D lists
x = exposure_times  # Independent variable
y = mean_values  # Dependent variable

# Reshape x to a 2D array because scikit-learn expects a 2D array for features
x_reshaped = np.array(x).reshape(-1, 1)
y = np.array(y)

# Create and fit the model
model = LinearRegression(fit_intercept=True)
model.fit(x_reshaped, y)

# Predict the values using the model
y_pred = model.predict(x_reshaped)

# Calculate residuals
residuals = y - y_pred

# Calculate Mean Squared Error
rmse = np.round(np.sqrt(mean_squared_error(y, y_pred)),3)

# Calculate R-squared
r_squared = np.round(r2_score(y, y_pred),2)

# Print the results
print(f"Coefficients: {model.coef_}")
print(f"Intercept: {model.intercept_}")
print(f"Root mean Squared Error: {rmse}")
print(f"R-squared: {r_squared}")

