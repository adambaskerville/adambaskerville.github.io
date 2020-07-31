#include <iostream>
#include <math.h>

// Define the function to be minimised, x**2 + y**2 + 3
double FuncValue (double x, double y)
{    
    return pow(x,2) + pow(y,2) + 3;
}

int main()
{	
    // Declare the data types of the input values x and y
    double x;
    double y;

    std::printf("This C++ program will be called from Python continuously\n");
    do 
    {   
        // This reads in the command line input and assigns to x and y
        std::cin >> x; 
        std::cin >> y; 
        
        // Print function value to the terminal
        std::printf("%.16f \n", FuncValue(x, y));

        // Clear cin
        std::cin.clear();
    } while (1 < 2); // Keep do loop running indefinitely

//while (std::cin.get() != 'Kill')
    return 0;
}