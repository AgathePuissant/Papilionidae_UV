
import cv2
import matplotlib.pyplot as plt
import numpy as np



def calculate_mean_rgb(image, positions):
    """
    Calculate mean RGB values for patches specified by positions.
    
    Args:
    image (numpy array): The image array.
    positions (list): List of tuples specifying the x, y, width, height of patches.
    
    Returns:
    list: A list of mean RGB values for each patch.
    """
    mean_rgb_values = []
    
    for (x, y, w, h) in positions:
        # Extract the patch
        patch = image[y:y+h, x:x+w]
        
        # Calculate mean RGB values (mean across the height and width dimensions)
        mean_b = np.mean(patch[:, :, 0])
        mean_g = np.mean(patch[:, :, 1])
        mean_r = np.mean(patch[:, :, 2])
        
        mean_rgb_values.append((mean_r, mean_g, mean_b))
    
    return mean_rgb_values

#Here put the colorchecker image OR the linearized image to check linearization
image_path = "path/to/colorchecker/image"
image = cv2.imread(image_path)
image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

# Define the positions of color checker patches in the NEF image
color_checker_positions = [
    (1875, 2275, 150, 150), 
    (1875, 1890, 150, 150), 
    (1875, 1560, 150, 150),
    (1875, 1200, 150, 150), 
    (1875, 860, 150, 150), 
    (1875, 500, 150, 150)
]
# Define the corresponding colors in RGB
reflectances = np.array([
    91.5736394356503,
    59.4147769174464,
    38.3956994557667,
    19.3768145989381,
    10.1723971251945,
    3.21941237074612
])


mean_rgb = calculate_mean_rgb(image, color_checker_positions)


plt.scatter(np.mean(mean_rgb,axis=1), reflectances)
