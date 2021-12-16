The file AVPR_A1_Klimack_Jason.m is the main file, which allows you 
to run a script, selecting a test image, and compute the nipple 
locations. 

The file nipple_klimack.m contains the function nipple_klimack, which 
takes a grayscale image and a boolean as input, and returns the same 
image, along with the x and y coordinates of the nipple locations. 
The input boolean, display, is used to toggle whether or not the result
of the nipple locations will be displayed on the screen as the 2 points
overlaid on the original image. 

nipple_klimack is called from within AVPR_A1_Klimack_Jason.m. 

There is an index value, idx, which allows control over which 
image to test. 

The PDF, AVPR_A1_Klimack.pdf is my report. 