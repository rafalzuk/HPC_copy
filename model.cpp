#include "MODEL.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

Model::Model(int argc, char* argv[], int rank, int size)
{
	ValidateParameters(argc, argv, rank, size); 
}

Model::Model()
{
}

Model::~Model()
{
}

void Model::ValidateParameters(int argc, char* argv[], int rank, int size)
{
	wrongParam = true; // parameters are apriori rejected unless all verified
	
	if(argc > 13)  // check if number of arguments is ok, only print message once
	{
		if(rank == 0)
		{
		cout << endl;
		cout << " Too many input arguments!" << endl;;
		RemindArguments();
		}
	} 
	else if(argc < 11) 
	{
		if(rank == 0)
		{
		cout << endl;
		cout << " Too few input arguments!" << endl;;
		RemindArguments();
		}
	} 
	else if(argc == 12) 
	{
		if(rank == 0)
		{
		cout << endl;
		cout << " Incorrect number of arguments!" << endl;;
		RemindArguments();
		}
	} 
	else
	{
		ss.str(argv[1]);
		if(ss >> ax) 
		{
			ss.clear();
			ss.str(argv[2]);
			if(ss >> ay) 
			{
				ss.clear();
				ss.str(argv[3]);
				if(ss >> b) 
				{
					ss.clear();
					ss.str(argv[4]);
					if(ss >> c) 
					{
						ss.clear();
						ss.str(argv[5]);
						if(ss >> T && T > 0) 
						{
							ss.clear();
							ss.str(argv[6]);
							if((ss >> inputTest && inputTest == ceil(inputTest)) && inputTest > 1) 
							{	
								Nt = (int)inputTest;
								ss.clear();
								ss.str(argv[7]);
								if(ss >> Lx && Lx > 0) 
								{
									ss.clear();
									ss.str(argv[8]);
									if((ss >> inputTest && inputTest == ceil(inputTest)) && inputTest > 1) 
									{
										Nx = (int)inputTest;
										ss.clear();
										ss.str(argv[9]);
										if(ss >> Ly && Ly > 0) 
										{
											ss.clear();
											ss.str(argv[10]);
											if((ss >> inputTest && inputTest == ceil(inputTest)) && inputTest > 1) 
											{
												
												Ny = (int)inputTest;
												if(argc == 13) // check whether the input contains extra data used in MPI
												{
													ss.clear();
													ss.str(argv[11]);
													if((ss >> inputTest && inputTest == ceil(inputTest)) && inputTest > 0)
													{
														Px = (int)inputTest;
														ss.clear();
														ss.str(argv[12]);
														if((ss >> inputTest && inputTest == ceil(inputTest)) && inputTest > 0)
														{
															ss.clear();
															Py = (int)inputTest;
															dx = Lx / (Nx - 1);
															dy = Ly / (Ny - 1);
															dt = T / (Nt - 1);
															if(Px <= Nx && Py <= Ny)
															{
																wrongParam = false;
															}
															
														}
													}
												}
												else // 11 parameters input i.e. we want single process
												{
													ss.clear();
													Px = 1;
													Py = 1;
													dx = Lx / (Nx - 1);
													dy = Ly / (Ny - 1);
													dt = T / (Nt - 1);
													wrongParam = false;
												}
											}
										}					
									}
								}
							}					
						}
					}
				}
			}
		}
	
		// At last check is the processes' number is ok
		if(!wrongParam && Px*Py != size)
		{
			if(rank == 0)
			{
				cout << endl;
				cout << " Number of processes specified in '-np' does not match the arguments input!" << endl;	
				cout << " To run the program with your own parameters use: mpiexec -np X YOUR ARGUMENTS" << endl;
				cout << " X should be equal to product of 11) * 12) arguments, or 1 if they do not exist!" << endl;
				RemindArguments();
			}
			wrongParam = true;
		}
		else if(wrongParam && rank == 0)
		{	
			cout << endl;
			cout << " One or more parameters incorrect!" << endl;
			RemindArguments();
		}
		else if(!wrongParam && rank == 0)
		{
			PrintParameters();
		}
	}
}

void Model::RemindArguments()
{
	cout << " You should provide 10 or 12 parameters of the following type and order:" << endl;
	cout << " 1) Real number for coefficient ax" << endl;
	cout << " 2) Real number for coefficient ay" << endl;
	cout << " 3) Real number for coefficient b" << endl;
	cout << " 4) Real number for coefficient c" << endl;
	cout << " 5) Real Positive number for max simulation time." << endl;
	cout << " 6) Integer > 1 for number of simulation time points." << endl;
	cout << " 7) Real Positive number for domain x width." << endl;
	cout << " 8) Integer > 1 for number of domain points along x dimension." << endl;
	cout << " 9) Real Positive number for domain y width." << endl;
	cout << "10) Integer > 1 for number of domain points along y dimension." << endl;
	cout << "	  ***INPUTS BELOW USED FOR MPI (OPTIONAL)***" << endl;
	cout << "11) Real Positive Integer for number of domain parts along x dimension." << endl;
	cout << "12) Real Positive Integer for number of domain parts along y dimension." << endl;
	cout << endl;
}

void Model::PrintParameters()
{
	cout << endl;
	cout << "Parameters accepted for a " << Px*Py << "-process computation."<< endl;
	cout << "The parameters of this model are: " << endl;
	cout << "ax    = " << ax << endl;
	cout << "ay    = " << ay << endl;
	cout << "b     = " << b << endl;
	cout << "c     = " << c << endl;
	cout << "T max = " << T << endl;
	cout << "Nt    = " << Nt << endl;
	cout << "Lx    = " << Lx << endl;
	cout << "Nx    = " << Nx << endl;
	cout << "Ly    = " << Ly << endl;
	cout << "Ny    = " << Ny << endl;
	cout << "dt    = " << dt << endl;
	cout << "dx    = " << dx << endl;
	cout << "dy    = " << dy << endl;
	if((Py > 1) || (Px > 1))
	{
	cout << "   MPI:" << endl;
	cout << "Px    = " << Px << endl;
	cout << "Py    = " << Py << endl;
	}
}

bool Model::IsValid()
{
	return !wrongParam;
}