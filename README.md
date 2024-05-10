#  Introduce
This is the code for the electromagnetic side of geophysical transients
The code consists of two main programs and five subroutines, namely: (1) numerical simulation main program, (2) measured data reading main program, (3) one-dimensional (uniform) media numerical simulation subroutine, (4) calculation of apparent resistance subroutine, (5) calculation of background resistivity subroutine, (6) correction of resistivity coefficient subroutine, and (7) time-depth conversion subroutine.
#  (1) numerical simulation main program
This main program is a numerical simulation script that calls the partial functions in the (3) to (7) subroutines to perform rapid imaging of homogeneous half-space or one-dimensional layered media. The script is divided into five parts, as outlined below:
1. The first part: Setting independent variables and model parameters. Manually set the geometry, physical property parameters and observation time window of the forward model.
2. The second part: one-dimensional (uniform) model forward modeling. Based on the parameters in the first part, the forward modeling mode (one-dimensional or uniform) is selected to get the electromagnetic field response.
3. The third part: Calculating the total apparent resistivity. Based on the second part, the electric field component EX and magnetic field component VZ are calculated in linear domain and logarithmic domain respectively.
4. The fourth part: resistivity correction. The whole period apparent resistivity obtained in the third part is corrected to eliminate geometric factors.
5. The fifth part: time-depth conversion. Based on the fourth part, the corresponding relation between resistivity and time is transformed into the corresponding relation between resistivity and depth.
#  (2) The measured data read main program
This main program is a script for fast imaging of the measured data. It will call the partial function in the subprogram (3) to (7) to complete the fast imaging of the measured data. The script is divided into four parts, as outlined below:
1. Part 1: Read data from txt file. Read in the measured data from the outside to prepare for subsequent imaging;
2. Part Two: Calculating the total apparent resistivity. Based on the second part, the electric field component EX and magnetic field component VZ are calculated in linear domain and logarithmic domain respectively.
3. The third part: resistivity correction. The whole period apparent resistivity obtained in the third part is corrected to eliminate geometric factors.
4. The fourth part: time-depth conversion. Based on the fourth part, the corresponding relation between resistivity and time is transformed into the corresponding relation between resistivity and depth.
#  (3) one-dimensional (uniform) media numerical simulation subroutine
This subroutine is a program for numerical simulation of one-dimensional or uniform media, in which only one subfunction (forwardmodeling) is used for numerical simulation, the input parameters are: geometric parameters, physical parameters, forwardmodeling model, and the output result is: electric field component EX and magnetic field component VZ. Note: The forwardmodeling requires the filter coefficient, and the filter coefficient file is named "filter coefficient". The reading position of the filter coefficient in the forwardmodeling subfunction should be reasonably modified in actual use.
#  (4) Calculating apparent resistance subroutine
This subroutine is used to calculate the lifetime apparent resistivity of the program, its input parameters are: field source and measuring point parameters, electromagnetic field response, output parameters are: linear domain and logarithmic and lifetime apparent resistivity. The program is composed of two parts: obtaining the peak time and calculating the resistivity, including nor_diff, log_diff, rou_ex and rou_vz. The following is a brief description of the subfunctions:
1.nor_diff: This subfunction is based on the linear domain and takes the peak time of the electromagnetic field components EX and VZ. ;
2.log_diff: This subfunction is based on the logarithmic domain to obtain the peak time of the electromagnetic field components EX and VZ;
3.rou_ex: This subfunction will call 1 and 2 above to calculate the total apparent resistivity of the electric field component EX;
4.rou_vz: This subfunction calls 1 and 2 above to calculate the lifetime apparent resistivity of the magnetic field component VZ;
#  (5) Calculating background resistivity subroutine
This subroutine is used to calculate the background resistivity of the program, its input parameters: field source and measuring point parameters, total apparent resistivity, output parameters: background resistivity value. The program consists of two parts: obtaining the background resistivity of numerical simulation and obtaining the background resistivity of measured data, including background_resistivity and background_resistivity_1 subfunctions. The following is a brief description of the subfunctions:
1. background_resistivity: This subfunction is used to calculate the background resistivity of the numerical simulation;
2. background_resistivity_1: This subfunction is used to calculate the background resistivity of the measured data;
#  (6) Correcting the resistivity coefficient subroutine
This subroutine is used to calculate the apparent resistivity correction coefficient of the program, its input parameters: field source and measuring point parameters, background resistivity and other parameters, output parameters: resistivity correction coefficient. The program consists of L_halfspace_T, rou_correction_1, rou_correction_2 and rou_correction_3. The following is a brief description of the subfunctions:
1. L_halfspace_T: This subfunction is used to calculate the electromagnetic field response of a uniform half-space long conductor source;
2. rou_correction_1: This subfunction is based on the uniform half-space forward program to calculate the correction coefficient;
3. rou_correction_2: This subfunction is based on one-dimensional media forward programming to calculate the correction coefficient;
4. rou_correction_3: This subfunction is to calculate the correction coefficient according to the measured data;
#  (7) Time depth conversion subroutine
This subroutine is used to carry out the time-depth transformation of the program, its input parameters: apparent resistivity parameters, observation time parameters, etc., output parameters: depth. The program contains two subfunctions T_H_1 and T_H_2. The following is a brief description of the subfunctions:
1. T_H_1: This subfunction is a time-depth transformation of the resistivity time relationship of the numerical simulation;
2. T_H_2: This subfunction is a time-depth transformation of the resistivity time relation of the measured data;
#  Requirements
Here are the packages in python that the program needs:
scipy.interpolate
numpy
pandas
scipy.constants
scipy.interpolate.CubicSpline
matplotlib.pyplot
math
time

The following is the filter coefficient file required for numerical simulation:
The name of the filtering coefficient is: Filtering coefficient.
