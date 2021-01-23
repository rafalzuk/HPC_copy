#include <chrono>
#include "BURGERS.h"
#include <ctime>

int main(int argc, char* argv[]) 
{
	int rank , size;
	int retval_rank, retval_size;
	int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) 
	{
        cout << " Failed to initialise MPI" << endl;
    }
	
	retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM)
	{
		cout << "Communicator invalid !" << endl;
		return 1;
	}

	Model model_1(argc, argv, rank, size); // create new model and check parameters
	
	if(model_1.IsValid()) // checke if inputs were fine
	{	
		
		Burgers burgers_1(model_1);
		burgers_1.SpaceVectors(rank, size);
		burgers_1.VelocityVectors(size);
		
		typedef std::chrono::high_resolution_clock hrc;
		typedef std::chrono::milliseconds ms;
		hrc::time_point start = hrc::now();
		
			burgers_1.Integrate(rank, size);
			
		hrc::time_point end = hrc::now();
		
		if(rank == 0)
		cout << "Integration complete" << endl;
		
		burgers_1.CalcEnergy(rank, size);
		burgers_1.WriteToFile(rank, size);
		
		
		cout << "Rank " << rank << " . Integration took: " << std::chrono::duration_cast<ms>(end - start).count() << " ms" << endl;
	}
	
	MPI_Finalize();
    return 0;
}
