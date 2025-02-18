# Low Filter Order Digital Equalization for Mobile Device Earphones

## Dependencies

[Discrete Frechet Distance Function](https://www.mathworks.com/matlabcentral/fileexchange/31922-discrete-frechet-distance)

[Balazs Bank Functions](https://home.mit.bme.hu/~bank/parfilt/#matlabcode)

From the first link download the function and place it inside the project folder.
From the second link download the `freqpoles`, `parfilt` and `parfiltidfr` functions and place them inside the project folder.

The `tfplots` function was originally developed by Balazs Bank and is included in the project folder since some modifications have been made to it.

## Description

This repository contains a demo version of the code used for the paper "Low Filter Order Digital Equalization for Mobile Device Earphones" by G. Kamaris, P. Zachos and J. Mourjopoulos.

You can run this using the `main.m` file and substitute the input and target functions with whatever you want. Apart from headphone equalization, this code can also be used for system identification purposes, by setting the target as the system to be identified and the input as the Dirac function. It has been used for this purpose in a paper that has been accepted in the Journal of the Audio Engineering Society and awaiting publication called 'Feedforward Headphone Active Noise Control Utilizing Auditory Masking'.

If you use any of this code in your work please cite the following paper:

Kamaris, G., Zachos, P., & Mourjopoulos, J. (2021). Low filter order digital equalization for mobile device earphones. Journal of the Audio Engineering Society, 69(5), 297-308.
