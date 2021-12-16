



% path to the test images
path = "test/test/";

% list of available test images
img_names = ["IR_0100.png", "IR_1367.png", "IR_2841.png", "IR_3635.png", "IR_3759.png", "IR_4089.png"];

% change idx to the index of the image you wish to test
idx = 1;
run_case(path + img_names(idx));


function run_case(filename)
    img = imread(filename);
    t1 = clock();
    % second argument is boolean to determine whether to display image
    overlay = nipple_klimack(img, 1);
    t2 = clock();
    t3 = t2 - t1;
    
    imgo = overlay{1};
    x = overlay{2};
    y = overlay{3};
	
end