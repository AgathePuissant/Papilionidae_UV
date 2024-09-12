import cv2
import numpy as np
import os
import csv

# Define the folder containing the images
image_folder = "path/to/raw/images"
output_csv = "path/to/output"

# Define the range for white color in HSV space
lower_white = np.array([0, 0, 200])
upper_white = np.array([180, 55, 255])

# Create a list to hold the detected coordinates
coordinates_list = []

# Iterate over all files in the folder
for filename in os.listdir(image_folder):
    if filename.lower().endswith(('.tif', '.png', '.jpg', '.jpeg', '.bmp', '.tiff')):
        image_path = os.path.join(image_folder, filename)
        
        # Load the image
        image = cv2.imread(image_path)
        
        if image is None:
            print(f"Failed to load image {filename}")
            continue
        
        # Convert the image to HSV color space
        hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
        
        # Create a mask for the white color
        mask = cv2.inRange(hsv, lower_white, upper_white)
        
        # Apply morphological operations to clean up the mask
        kernel = np.ones((5, 5), np.uint8)
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
        
        # Find contours in the mask
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        # Process the contours
        for contour in contours:
            area = cv2.contourArea(contour)
            if area > 500:  # Adjust the area threshold as needed
                # Approximate the contour to a circle
                (x, y), radius = cv2.minEnclosingCircle(contour)
                center = (int(x), int(y))
                radius = int(radius)
                
                # Save the detected coordinates
                coordinates_list.append({'filename': filename, 'x': center[0], 'y': center[1]})
                
                # Optionally, draw the circle in the original image for verification
                cv2.circle(image, center, radius, (0, 255, 0), 2)
                cv2.rectangle(image, (center[0] - 5, center[1] - 5), (center[0] + 5, center[1] + 5), (0, 128, 255), -1)
                
                # Save or display the image with the detected circle (for debugging)
                # cv2.imwrite(os.path.join('output_folder', filename), image)
                # cv2.imshow("Detected Circle", image)
                # cv2.waitKey(0)
                
# Write the coordinates to a CSV file
with open(output_csv, mode='w', newline='') as csvfile:
    fieldnames = ['filename', 'x', 'y']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    for row in coordinates_list:
        writer.writerow(row)

print(f"Coordinates saved to {output_csv}")
