function[]=problem_set1()





%%%%%%%%%%%%%%%problem 1: load and plot %%%%%%%%%%%%%%
%%%%%COMMANDS: load, size, plot

  %load in a file, data_prob2.dat. It is xy data, column 1 is x, column 2 is y


%what is the size of this?

%how many rows?

%how many columns?

%rename column1 as xvar

%rename column2 as yvar

%multiply column 2 by 2, name as y2.

%plot xvar vs the new variable, y2.




%%%%%%%%problem 2, function input output %%%%%%%%%
%%%%COMMANDS: plot, xlabel, ylabel

%create a new script with the first line: 
  %function[xValues, yValues]=calc_gauss(mu, sigma)

%define xmin as mu-sigma*4

%define xmax=mu+sigma*4

  %define nPts = 100;
  %define dx=(xmax-xmin)/100

%create an array of xValues from xmin to xmax, at increments of dx

%define the gaussian f=1/sqrt(2*pi*sigma^2)*exp(-1/(2*sigma^2)*(x-mu)^2)

%return the xValues and yValues that were calculated in this function

  %NOTE: functions in MATLAB are more like c++ than python. If you call a function, it will not
%have access to anything outside of that function unless you pass it as input. The code that calls that function
%will not have access to anything computed by that function unless it is returned as output.

  %%%%%%%%%% problem 3: call your function from another program %%%
%%%%%%plot and edit the plot.
  %%% To call a function that is defined externally, there are a few options. If the function is defined in a distinct file, and is in the same directory, you can call it. Alternatively, in MATLAB you can define the function at the bottom of the 'main' file that is calling it.  
%open a new matlab script
  %You can also create libraries that contain function definitions, and then include or in python, import them. That is a bit more trouble, but can be worth it if you want to re-use the functions in many other programs. 

  %Create a for loop from 1:10

					   %define width= i*0.1, so it increases as you iterate

%Generate a gaussian using the calc_gauss function you wrote in problem 3. 

%Use the same mean value mu for each iteration

%Feed in the new width variable for sigma, so the variance keeps growing with each iteration of the for loop

%plot x vs f for each value of sigma.

					   %Make sure that you "hold" the figure before you enter the for loop, so that each curve gets added on to the same figure.

%It should look like a spectrum of normal distributions that keep getting wider and lower. 

%label xaxis x

%label yaxis p(x)

%change the fontsize to 40

					   %After adding all the curves,
%save the figure using the saveas command to a file named 'myGaussians.eps'




%%%%%%%%%%%%problem 4, surface & contour plots %%%%%%%%%%%%
%%%%COMMANDS: meshgrid, surf, contour

%create a meshgrid of x values from -5 to +5, and y values from -10 to +10

%define a function on this grid Z=X^2*cos(X)-Y^2;

%create a surface plot of Z vs X and Y

% hold this plot

%add a contour plot of Z vs X and Y



%%%%%%%%problem 5, calculate means and errors, plot%%%%%%%
%%%%%COMMANDS: load, size, mean, std, errorplot 

%load datapts3.dat file

%what is the size of the array?

%calculate the mean of the data pts

%calculate the standard deviation of the data pts

%calculate the standard error of the mean (SEM) for the first 10 data points.
%SEM is the sample standard deviation divided by the square root of the
%sample size.

%calculate the mean of the first 10 data points

%create an array that stores the number of data points you just sampled. (i.e.
%npts(1)=10;

%create an array to store the SEM you just calculated (i.e. semvec(1)=..)
%and another for the mean mymeanvec(1)=..)


%calculate the standard error of the mean for the first 100 data points.

%store these as the next entries in your datapoints array and SEM array

%repeat for the first 200 and 1000 datapoints

%make an errorplot of the mean as a function of number of datapoints, with
%the SEM used as the size of the bars



   %%%%%%%%%%%%problem 6. 
 %  Download Visual Molecular Dynamics (VMD) for free from:  www.ks.uiuc.edu)
   %%%load in a file (new molecule). 
 %%load in the clath7_col.psf. VMD should recognize it as a PSF file.
   %%Then load in the trajectory: clathrin713frames.xyz
%%   Under graphics->Representations, change the coloring method to a different ColorID (not gray). Change the drawing method to licorice, with a bond radius of 1.7. Watch the movie, and Render (under File) an image of one of the frames.



   %%%%%%%%%%problem 7, physical biology of the cell, problem 2.6%%%
%%%%%%estimating organelle sizes%%%%%
