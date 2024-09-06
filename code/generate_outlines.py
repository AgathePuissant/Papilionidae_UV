# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import tqdm


# Function to pad image with white pixels to target size
def pad_image(img, target_size):
    padded_img = np.ones((target_size[1], target_size[0], 3), dtype=np.uint8) * 255  # Create a white image
    h, w = img.shape[:2]
    y_off = (target_size[1] - h) // 2
    x_off = (target_size[0] - w) // 2
    padded_img[y_off:y_off+h, x_off:x_off+w, :] = img  # Place the original image in the center
    return padded_img


path_im=r"C:\Users\Agathe\Desktop\visible_pictures"
path_im_save=r"C:\Users\Agathe\Desktop\visible_pictures\visible_pictures_resized"
path_outlines = r"C:\Users\Agathe\Desktop\visible_pictures\outlines_visible"
listim = os.listdir(path_im)
listim = [x for x in listim if x[-4:]==".jpg"]

# Load your images and determine the maximum dimensions
max_width = 0
max_height = 0
images = []
for image_path in listim:
    img = cv2.imread(path_im+"/"+image_path)
    max_width = max(max_width, img.shape[1])
    max_height = max(max_height, img.shape[0])
    images.append(img)

# Pad images to the maximum dimensions
padded_images = []
for img in tqdm.tqdm(images):
    padded_img = pad_image(img, (max_width, max_height))
    padded_images.append(padded_img)

# Now all images in padded_images will have the same size, padded with white pixels


for i in tqdm.tqdm(range(len(padded_images))) :
    im= padded_images[i]
    
    im=cv2.resize(im, (im.shape[1]//4,im.shape[0]//4)) #Here you can resize the image for shorter computation time for the patternize code
    
    gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    
    grayblur = cv2.blur(gray, [10,10])
    
    ret, thresh = cv2.threshold(grayblur, 250,255, cv2.THRESH_BINARY_INV)
    
    thresh= cv2.erode(thresh, np.ones([15,15]))

    # Find contours
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    contours = sorted(contours, key=cv2.contourArea, reverse=True)[0]
    
    xy=contours.squeeze()
    np.savetxt(path_outlines+"/"+listim[i][:-4]+".txt", xy, fmt='% 4d')
    
    cv2.imwrite(path_im_save+"/"+listim[i], im)
