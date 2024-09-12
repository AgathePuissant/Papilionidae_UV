import cv2
import numpy as np
import os
from glob import glob
import csv

# Affine regression parameters
a = 2.35558407
b = 5.4339752588951455

# Define paths
input_folder = "path/to/raw/images"  # Replace with your input folder path
output_folder = "path/to/output/folder"  # Replace with your output folder path
coordinates_csv = "path_to_spectralon/coordinates/for/each/image"  # Path to your CSV file with coordinates

# Ensure the output directory exists
os.makedirs(output_folder, exist_ok=True)

# Read the coordinates from the CSV file
coordinates_dict = {}

with open(coordinates_csv, mode='r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        filename = row['filename']
        x = int(row['x'])
        y = int(row['y'])
        coordinates_dict[filename] = (x, y)

# Function to calibrate an image using the white standard area
def calibrate_image_using_white_standard(image, a, b, white_standard_value):
    """
    Calibrates the image using the white standard value and the affine relationship.
    
    Args:
        image (numpy.ndarray): Grayscale image to be calibrated.
        a (float): Slope of the affine regression (P = aR + b).
        b (float): Intercept of the affine regression (P = aR + b).
        white_standard_value (float): Pixel value of the white standard region.
    
    Returns:
        numpy.ndarray: Calibrated image with reflectance values.
    """
    # Convert the image to float for accurate calibration
    image_float = image.astype(np.float32)
    
    # Calculate the reflectance of the white standard
    R_white = (white_standard_value - b) / a
    
    
    if R_white == 0:
        raise ValueError("Reflectance of the white standard is zero, which will cause division by zero.")
    
    # Apply the inverse affine transformation to the entire image
    reflectance = (image_float - b) / a
    
    # Normalize the reflectance by the reflectance of the white standard
    reflectance /= R_white
    
    # Ensure reflectance is within the valid range [0, 1]
    reflectance = np.clip(reflectance, 0, 1)
    
    return reflectance

# Process each TIFF image in the input folder
for image_path in glob(os.path.join(input_folder, '*.tif')):
    filename = os.path.basename(image_path)
    
    # Check if coordinates for this image are available
    if filename not in coordinates_dict:
        print(f"No coordinates found for {filename}. Skipping.")
        continue

    x, y = coordinates_dict[filename]
    w, h = 50, 50  # Width and height for the white standard area

    # Load the grayscale image
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    
    if image is None:
        print(f"Failed to load image {image_path}")
        continue

    # Ensure the white standard area is within the image bounds
    if (x + w > image.shape[1]) or (y + h > image.shape[0]):
        print(f"White standard area for {filename} is out of image bounds. Skipping.")
        continue
    
    # Extract the white standard area
    white_standard = image[y:y+h, x:x+w]

    # Calculate the average intensity of the white standard
    avg_intensity = white_standard.mean()

    if avg_intensity == 0:
        print(f"White standard in {filename} has zero intensity. Skipping image.")
        continue

    # Calibrate the image using the white standard value
    try:
        calibrated_image = calibrate_image_using_white_standard(image, a, b, avg_intensity)

        # Convert reflectance to pixel values in the range [0, 255]
        calibrated_image = (calibrated_image * 255).astype(np.uint8)

        # Generate the output path
        output_image_path = os.path.join(output_folder, filename)

        # Save the calibrated image
        cv2.imwrite(output_image_path, calibrated_image)

        print(f"Processed and saved {output_image_path}")
    except ValueError as e:
        print(f"Error processing {filename}: {e}")

print("All images processed.")
