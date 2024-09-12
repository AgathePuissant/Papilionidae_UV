import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from numpy.polynomial.polynomial import Polynomial

# Function to extract RGB values from color checker patches in an image
def extract_color_checker_rgb(image_path, color_checker_positions):
    image = cv2.cvtColor(cv2.imread(image_path, cv2.IMREAD_UNCHANGED), cv2.COLOR_BGR2RGB) / 65535.0
    rgb_values = []
    for (x, y, w, h) in color_checker_positions:
        patch = image[y:y+h, x:x+w]
        avg_rgb = np.mean(patch, axis=(0, 1))
        rgb_values.append(avg_rgb)
    return np.array(rgb_values)

# Function to fit a polynomial for each RGB channel
def fit_polynomial_model(captured_rgb, reference_rgb, degree=3):
    models = []
    for i in range(3):  # Separate model for R, G, and B channels
        coeffs = np.polynomial.polynomial.Polynomial.fit(captured_rgb[:, i], reference_rgb[:, i], degree).convert().coef
        models.append(coeffs)
    return models

# Function to apply the polynomial model to an image
def apply_polynomial_model(image, models):
    corrected_image = np.zeros_like(image)
    for i in range(3):  # Apply model to R, G, and B channels independently
        poly = Polynomial(models[i])
        corrected_image[..., i] = poly(image[..., i])
    return np.clip(corrected_image, 0, 1)

# Define the positions of color checker patches in the image
color_checker_positions = [
    (1875, 2275, 150, 150), 
    (1875, 1890, 150, 150), 
    (1875, 1560, 150, 150),
    (1875, 1200, 150, 150), 
    (1875, 860, 150, 150), 
    (1875, 500, 150, 150)
]

# Reference RGB values for the gray patches (normalized)
reference_rgb_values = np.array([
    [91.5736394356503, 91.5736394356503, 91.5736394356503],
    [59.4147769174464, 59.4147769174464, 59.4147769174464],
    [38.3956994557667, 38.3956994557667, 38.3956994557667],
    [19.3768145989381, 19.3768145989381, 19.3768145989381],
    [10.1723971251945, 10.1723971251945, 10.1723971251945],
    [3.21941237074612, 3.21941237074612, 3.21941237074612]
]) / 100

# Root directory containing the image files
root_dir = "path/to/images"

# List of image files to process
nef_files = glob(os.path.join(root_dir, "*.tif"))

output_folder = "path/to/output/folder"
os.makedirs(output_folder, exist_ok=True)

# Process each image
for nef_file in nef_files:
    print(f"Processing: {nef_file}")
    
    # Extract RGB values from the color checker
    captured_rgb_values = extract_color_checker_rgb(nef_file, color_checker_positions)
    
    # Fit polynomial models for each channel
    models = fit_polynomial_model(captured_rgb_values, reference_rgb_values)
    
    # Read the image file
    image = cv2.imread(nef_file, cv2.IMREAD_UNCHANGED)
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB) / 65535.0
    
    # Apply the polynomial model for color correction
    corrected_image = apply_polynomial_model(image, models)
    corrected_image_srgb = corrected_image * 65535.0
    
    # Convert back to uint8 for saving
    corrected_image_uint8 = corrected_image_srgb.astype(np.uint8)
    
    # Save the corrected image
    output_path = os.path.join(output_folder, os.path.basename(nef_file))
    cv2.imwrite(output_path, corrected_image_uint8)
    print(f"Saved corrected image: {output_path}")

print("Color correction applied to all images.")