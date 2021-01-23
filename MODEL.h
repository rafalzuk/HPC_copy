#ifndef CLASS_MODEL
#define CLASS_MODEL
#include <mpi.h>
#include <cmath>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <limits>
#include <iostream>
#include <string>
#include <array>
#include <sstream>
#include <istream>
#include <fstream>
using namespace std;

class Model
{
	public:
		Model(int argc, char* argv[], int rank, int size);
		Model();
        ~Model();
		
		void RemindArguments();
        void PrintParameters();
        bool IsValid();

        // Getters
		bool   IsVerbose() const { return verbose; }
        bool   IsHelp()    const { return help; }
        double GetX0()     const { return x0; }
        double GetY0()     const { return y0; }
        double GetLx()     const { return Lx; }
        double GetLy()     const { return Ly; }
        double GetT()      const { return T; }
        int    GetNx()     const { return Nx; }
        int    GetNy()     const { return Ny; }
        int    GetNt()     const { return Nt; }
        double GetDx()     const { return dx; }
        double GetDy()     const { return dy; }
        double GetDt()     const { return dt; }
        double GetAx()     const { return ax; }
        double GetAy()     const { return ay; }
        double GetB()      const { return b; }
        double GetC()      const { return c; }
		int GetPx()		   const {return Px; }
		int GetPy()		   const {return Py; }

		
	private:
        void ParseParameters(int argc, char* argv[]);
        void ValidateParameters(int argc, char* argv[], int rank, int size);

		bool wrongParam;
	    bool verbose;
        bool help;

	    // Numerics
	    double x0 = 0;
	    double y0 = 0;
	    double Lx;
	    double Ly;
	    double T;
	    int    Nx;
	    int    Ny;
	    int    Nt;
	    double dx;
	    double dy;
	    double dt;
		double inputTest;
		
		int Px;
		int Py;
		//////////
		
	    // Physics
	    double ax;
	    double ay;
	    double b;
	    double c;
		
        // Add any additional parameters here...
		std::string s1;
		std::stringstream ss;
};

#endif