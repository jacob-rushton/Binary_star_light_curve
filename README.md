# Binary_star_light_curve
A collection of files used to analyze optical images taken by the Plaskett Telescope of an eclipsing binary star system referred to in literature as 1SWASP J074658.62+224448.5

This is a very complicated code involving a some astronomical imaging jargon, so let me try and and explain what the code does step by step.

## 1: Calibration
The first step involves calibrating each optical image, or "science" frame;  that involves trimming any non-information out of the image ad then using the flat and dark frames to calibrate each image.  Combining all dark frames gives a "raw master dark", which is then subtracted from each flat frame.  These subtracted flat frames are then combined into a "master flat".  Then the raw master dark pixel values are subtracted from each science frame and then the pixels of each resulting science frame are then divided by the corresponding pixel of the master flat. The frames are now calibrated.

## 2: Alignment
The imaging took place over two days, and each day the imaging went on for a period of hours.  To correct for the movement of the Earth, the stars in the images must be aligned so that the pixels containing the stars of interest are maintained
