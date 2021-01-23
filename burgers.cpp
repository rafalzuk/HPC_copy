#include "BURGERS.h"
#include "cblas.h"
//#include <array>
//#include <numeric>

Burgers::Burgers(const Model& model)
{

	/*Import some parameters, to avoid calling getters 
	of Model class every time these are needed, thus boost performance */
	Dx = model.GetDx();
	Dy = model.GetDy();
	Dt = model.GetDt();
	
	// Physical Constants
	ax = model.GetAx();
	ay = model.GetAy();
	b  = model.GetB();
	c  = model.GetC();
	
	/* Values determining size of arrays
	and time integration points */
	sizeX = model.GetNx();
	sizeY = model.GetNy();
	sizeXY = sizeX*sizeY;
	Px = model.GetPx();
	Py = model.GetPy();
	Nt = model.GetNt();
	
	/* Halves of the domain size */
	Lx_half = model.GetLx() / 2;
	Ly_half = model.GetLy() / 2;
	
}

// Defalut constructor
Burgers::Burgers()
{
	cout << "Please provide a model as an argument! " << endl;
}

// Destructor
Burgers::~Burgers()
{
	
}

void Burgers::SpaceVectors(int rank, int size)
{	
			// Build local X vectors for each process
			if(sizeX%Px == 0) // If divides evenly
			{
				local_sizeX = (sizeX/Px);
				const int LOCALSIZEX = local_sizeX;
				X_local = new double[LOCALSIZEX];
				for (int j = 0 ; j < local_sizeX ; j++ )
				{  
					X_local[j] = -Lx_half + j*Dx + (rank%Px)*(sizeX/Px)*Dx;
				}
			}
			else // Divides with remainder
			{
				if((rank+1)%Px != 0) // If is not the last process in the x dimension
				{
					local_sizeX = (int)(sizeX/Px);
					const int LOCALSIZEX = local_sizeX;
					X_local = new double[LOCALSIZEX];
					for (int j = 0 ; j < local_sizeX ; j++ )
					{  
						X_local[j] = -Lx_half + j*Dx + (rank%Px) *local_sizeX * Dx;
					}
				}
				else // If the last process in the x dimension
				{
					local_sizeX = sizeX - (int)(sizeX/Px)*(Px-1);
					const int LOCALSIZEX = local_sizeX;
					X_local = new double[LOCALSIZEX];
					for (int j = 0 ; j < local_sizeX ; j++ )
					{  
						X_local[j] = -Lx_half + j*Dx + (rank%Px) * (int)(sizeX/Px) * Dx;
					}
				}
			}
			
			
			// Similarly build local Y vector for each process
			if(sizeY%Py == 0)
			{
				local_sizeY = (sizeY/Py);
				const int LOCALSIZEY = local_sizeY;
				Y_local = new double[LOCALSIZEY];
				for (int j = 0 ; j < local_sizeY ; j++ )
				{  
					Y_local[j] = -Ly_half + j*Dy + (rank/Px)*local_sizeY*Dy;
				}
			}
			else
			{
				if( ((rank/Px)+1) != Py)
				{
					local_sizeY = (int)(sizeY/Py);
					const int LOCALSIZEY = local_sizeY;
					Y_local = new double[LOCALSIZEY];
					for (int j = 0 ; j < local_sizeY ; j++ )
					{  
						Y_local[j] = -Ly_half + j*Dy + (rank/Px)*local_sizeY * Dy;
					}
				}
				else
				{
					local_sizeY = sizeY - (int)(sizeY/Py)*(Py-1);
					const int LOCALSIZEY = local_sizeY;
					Y_local = new double[LOCALSIZEY];
					for (int j = 0 ; j < local_sizeY ; j++ )
					{  
						Y_local[j] = -Ly_half + j*Dy + (rank/Px)*(int)(sizeY/Py) * Dy;
					}
				}
			}
			
	
}

void Burgers::VelocityVectors(int size)
{	
	/* Conditions below determine the outermost points 
	 * (counting in absolute values i.e. from 1)
	 * which might be affected by initial condition
	 * i.e. become non-zero in unit central circle.
	 * This helps to reduce computational power required
	 * to implement initial veolcity field, as computer won't 
	 * need to compute radius everywhere in the domain.
	 */
	if(size == 1)
	{
	
		if (Lx_half > 1 && (sizeX%2 == 0))
		{
			max_xID = sizeX/2 + ceil(1/Dx);
			min_xID = sizeX/2 + 1 - ceil(1/Dx);
		}
		else if (Lx_half > 1 && (sizeX%2 != 0))
		{
			max_xID = (sizeX+1)/2 + ceil(1/Dx);
			min_xID = (sizeX+1)/2 - ceil(1/Dx);
		}
		else 
		{
			max_xID = sizeX;
			min_xID = 1;
		}
		
		if (Ly_half > 1 && (sizeY%2 == 0))
		{
			max_yID = sizeY/2 + ceil(1/Dy);
			min_yID = sizeY/2 + 1 - ceil(1/Dy);
		}
		else if (Ly_half > 1 && (sizeY%2 !=0))
		{
			max_yID = (sizeY+1)/2 + ceil(1/Dy);
			min_yID = (sizeY+1)/2 - ceil(1/Dy);
		}
		else
		{
			max_yID = sizeY;
			min_yID = 1;
		}
		
		const int SizeXY = sizeXY;
		U = new double[SizeXY];
		V = new double[SizeXY];
		U_next = new double[SizeXY];
		V_next = new double[SizeXY];
		
		// intialise velocity field by putting zeros
		for (int i = 0; i < sizeXY; i++)
		{
			U[i] = 0;
			V[i] = 0;
			U_next[i] = 0;
			V_next[i] = 0;
		}
		
		/* Now assign initial velocity field in the reduced domain
		 * noting that array indexes are counted with 0 and max/min_ID's
		 * were with 1 for convenience*/
		for(int j = (min_yID - 1); j < max_yID ; j++ )
		{ 
			for(int i = (min_xID - 1); i < max_xID ; i++ )
			{
				radius = pow( X_local[i]*X_local[i] + Y_local[j]*Y_local[j] , 0.5);

				if (radius <= 1.0)
				{
					// These are unit circle Velocities at time 0
					U[i+j*sizeX] = 2 * pow( (1-radius) , 4) * (4*radius + 1);
					V[i+j*sizeX] = 2 * pow( (1-radius) , 4) * (4*radius + 1);
					
				}
				
			}
			
		}
		
	}
	else // Now solution for multi-process computation
	{	
		const int localSize = local_sizeX*local_sizeY;
		U = new double[localSize];
		V = new double[localSize];
		U_next = new double[localSize];
		V_next = new double[localSize];
		
		
		for(int j = 0; j < local_sizeY ; j++ )
		{ 
			for(int i = 0; i < local_sizeX ; i++ )
			{
				U_next[i+j*local_sizeX] = 0;
				V_next[i+j*local_sizeX] = 0;
				
				radius = pow( X_local[i]*X_local[i] + Y_local[j]*Y_local[j] , 0.5);
				
				if (radius <= 1.0)
				{
					// These are unit circle Velocities at time 0
					U[i+j*local_sizeX] = 2 * pow( (1-radius) , 4) * (4*radius + 1);
					V[i+j*local_sizeX] = 2 * pow( (1-radius) , 4) * (4*radius + 1);
				}
				else
				{
					U[i+j*local_sizeX] = 0;
					V[i+j*local_sizeX] = 0;
				}
				
			}
			
		}
		
	}
	
}

void Burgers::Integrate(int rank, int size)
{	
	if(size == 1) // This single-process case exists to avoid if statements present in multi-process part
	{
		for(int t = 1; t < Nt; t++) // for every time point, after 0th which was defined above
		{
			for(int y = 1; y < (sizeY-1) ; y++) // for every y direction point, execpt for boundaries
			{
				for(int x = 1; x < (sizeX-1) ; x++) // for every x direction point, execpt for boundaries
				{
					// 'next' refers to the next time point
					U_next[x+y*sizeX] = Dt * (
									    c*( U[x+1+y*sizeX] - 2*U[x+y*sizeX] + U[x-1+y*sizeX] ) / (Dx*Dx) + 
										c*( U[x+(y+1)*sizeX] - 2*U[x+y*sizeX] + U[x+(y-1)*sizeX] ) / (Dy*Dy) -
										(ay + b*V[x+y*sizeX]) * (U[x+y*sizeX] - U[x+(y-1)*sizeX]) / Dy -
										(ax + b*U[x+y*sizeX]) * (U[x+y*sizeX] - U[x-1+y*sizeX]) / Dx 
											  ) + U[x+y*sizeX];
										
					V_next[x+y*sizeX] = Dt * ( 
										c*( V[x+1+y*sizeX] - 2*V[x+y*sizeX] + V[x-1+y*sizeX] ) / (Dx*Dx) + 
										c*( V[x+(y+1)*sizeX] - 2*V[x+y*sizeX] + V[x+(y-1)*sizeX] ) / (Dy*Dy) -
										(ay + b*V[x+y*sizeX]) * (V[x+y*sizeX] - V[x+(y-1)*sizeX]) / Dy -
										(ax + b*U[x+y*sizeX]) * (V[x+y*sizeX] - V[x-1+y*sizeX]) / Dx 
											  ) + V[x+y*sizeX];		
				}
			}
						
			// Prepare arrays for next time iteration
			for (int i = 0; i < sizeXY; i++)
			{
				U[i] = U_next[i];
				V[i] = V_next[i];
			}
			
		}
	
	} //Close single process procedure
	else // Now consider multi-process computation
	{
		// First check for boundary processes
		if(rank%Px == 0) // processes on the left edge of domain
		{
			if(rank == 0) //bottom left corner
			{
				// If the number of processses in a dimension is one, then do not send anything in this dimension
				// Also adjust local domain size accordingly for sending data in the perpendicular directoin
				tempYsize = (Py == 1) ? (local_sizeY-1) : local_sizeY;
				if(Px == 1)
				{
					tempXsize = (local_sizeX-1);
				}
				else
				{
					tempXsize = (local_sizeX);
					//Variables used for data sending and receiving
					U_right_send = new double[tempYsize-1];
					V_right_send = new double[tempYsize-1];
					U_right_import = new double[tempYsize-1];
					V_right_import = new double[tempYsize-1];
				}
				if(Py != 1)
				{
					//Variables used for data sending and receiving
					U_up_send = new double[tempXsize-1];
					V_up_send = new double[tempXsize-1];
					V_up_import = new double[tempXsize-1];
					U_up_import = new double[tempXsize-1];
				}
				
				for(int T = 1; T < Nt; T++)
				{
					for(int y = 1; y < tempYsize; y++)
					{
						for (int x = 1; x < tempXsize ; x++)
						{
							// Input which may change depending on process location in the domain
							if(x == 1 && y == 1 && Px != 1) // exchange data with right process
							{								
								// Create data packagae to be sent
								for(int i = 0; i < tempYsize-1; i++)
								{
									U_right_send[i] = U[(y+i)*local_sizeX+(local_sizeX-1)];
									V_right_send[i] = V[(y+i)*local_sizeX+(local_sizeX-1)];
								}
								
								// Send data
								MPI_Send(U_right_send, tempYsize-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
								MPI_Send(V_right_send, tempYsize-1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_right_import, tempYsize-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_right_import, tempYsize-1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							// If on the boundary, use received data
							if(x == (local_sizeX-1))
							{
									u_right = U_right_import[y-1];
									v_right = V_right_import[y-1];
							}
							else{
									u_right = U[x+1+y*local_sizeX];
									v_right = V[x+1+y*local_sizeX];
							}
							
							if(y == 1 && x == 1 && Py != 1) // exchange data with above process
							{								
								for(int i = 0; i < tempXsize-1; i++)
								{
									U_up_send[i] = U[(local_sizeY-1)*local_sizeX+i+x];
									V_up_send[i] = V[(local_sizeY-1)*local_sizeX+i+x];
								}
								
								MPI_Send(U_up_send, tempXsize-1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_up_send, tempXsize-1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_up_import, tempXsize-1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_up_import, tempXsize-1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								
							}
							if(y == (local_sizeY-1))
							{
									u_up = U_up_import[x-1];
									v_up = V_up_import[x-1];
							}
							else{
								u_up = U[x+(1+y)*local_sizeX];
								v_up = V[x+(1+y)*local_sizeX];
							}
							
							// Now integrate
							U_next[x+y*local_sizeX] = Dt * (
										c*( u_right - 2*U[x+y*local_sizeX] + U[x-1+y*local_sizeX] ) / (Dx*Dx) + 
										c*( u_up - 2*U[x+y*local_sizeX] + U[x+(y-1)*local_sizeX] ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - U[x+(y-1)*local_sizeX]) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - U[x-1+y*local_sizeX]) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( v_right - 2*V[x+y*local_sizeX] + V[x-1+y*local_sizeX] ) / (Dx*Dx) + 
										c*( v_up - 2*V[x+y*local_sizeX] + V[x+(y-1)*local_sizeX] ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] - V[x+(y-1)*local_sizeX]) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - V[x-1+y*local_sizeX]) / Dx 
											  ) + V[x+y*local_sizeX];	
	
						}
					}
					// Prepare vector for next time interation
					for (int i = 0; i < local_sizeX*local_sizeY; i++)
					{
						U[i] = U_next[i];						
						V[i] = V_next[i];
					}
					
				}
				
			} // Close bottom left corner
			else if(rank == Px*(Py-1)) // Process in the top left corner of the domain
			{
				if(Px == 1)
				{
					tempXsize = (local_sizeX-1);
				}
				else
				{		
					tempXsize = (local_sizeX);	
					U_right_send = new double[local_sizeY-1];
					V_right_send = new double[local_sizeY-1];
					U_right_import = new double[local_sizeY-1];
					V_right_import = new double[local_sizeY-1];
				}
				
				U_down_send = new double[tempXsize-1];
				V_down_send = new double[tempXsize-1];
				U_down_import = new double[tempXsize-1];
				V_down_import = new double[tempXsize-1];
				
				for(int T = 1; T < Nt; T++)
				{
					for (int y = 0; y < local_sizeY - 1 ; y++)
					{
						for(int x = 1; x < tempXsize; x++)
						{							
							if(x == 1 && y == 0 && Px != 1) // exchange data with right process
							{								
								for(int i = 0; i < local_sizeY-1; i++)
								{
									U_right_send[i] = U[(y+i)*local_sizeX+(local_sizeX-1)];
									V_right_send[i] = V[(y+i)*local_sizeX+(local_sizeX-1)];
								}
									
								MPI_Send(U_right_send, local_sizeY-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
								MPI_Send(V_right_send, local_sizeY-1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_right_import, local_sizeY-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_right_import, local_sizeY-1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							
							if(x == (local_sizeX-1))
							{
									u_right = U_right_import[y];
									v_right = V_right_import[y];
							}
							else{
									u_right = U[x+1+y*local_sizeX];
									v_right = V[x+1+y*local_sizeX];
							}
							
					
							if(y == 0 && x == 1) // exchange data with below process
							{
								for(int i = 0; i < tempXsize-1; i++)
								{
									U_down_send[i] = U[x+i];
									V_down_send[i] = V[x+i];
								}
								
								MPI_Send(U_down_send, tempXsize-1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_down_send, tempXsize-1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_down_import, tempXsize-1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_down_import, tempXsize-1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								
							}
							if(y == 0)
							{
									u_down = U_down_import[x-1];
									v_down = V_down_import[x-1];
							}
							else{
									u_down = U[x+(y-1)*local_sizeX];
									v_down = V[x+(y-1)*local_sizeX];
							}
							
							U_next[x+y*local_sizeX] = Dt * (
										c*( u_right - 2*U[x+y*local_sizeX] + U[x-1+y*local_sizeX] ) / (Dx*Dx) + 
										c*( U[x+(y+1)*local_sizeX] - 2*U[x+y*local_sizeX] + u_down ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - U[x-1+y*local_sizeX]) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( v_right - 2*V[x+y*local_sizeX] + V[x-1+y*local_sizeX] ) / (Dx*Dx) + 
										c*( V[x+(y+1)*local_sizeX] - 2*V[x+y*local_sizeX] + v_down) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - V[x-1+y*local_sizeX]) / Dx 
											  ) + V[x+y*local_sizeX];	
						}
					}
					
					for (int i = 0; i < local_sizeX*local_sizeY; i++)
					{
						U[i] = U_next[i];						
						V[i] = V_next[i];
					}
				}
			} // Close top left corner
			else // Remaining Processes in the left side of domain (w/o corner processes)
			{	
				if(Px == 1)
				{
					tempXsize = (local_sizeX-1);
				}
				else
				{		
					tempXsize = (local_sizeX);	
					U_right_send = new double[local_sizeY];
					V_right_send = new double[local_sizeY];
					U_right_import = new double[local_sizeY];
					V_right_import = new double[local_sizeY];
				}
				
				U_down_send = new double[tempXsize-1];
				V_down_send = new double[tempXsize-1];
				U_down_import = new double[tempXsize-1];
				V_down_import = new double[tempXsize-1];
				U_up_send = new double[tempXsize-1];
				V_up_send = new double[tempXsize-1];
				U_up_import = new double[tempXsize-1];
				V_up_import = new double[tempXsize-1];
				
				for(int T = 1; T < Nt; T++)
				{
					for (int y = 0; y < local_sizeY ; y++)
					{
						for(int x = 1; x < tempXsize; x++)
						{
							
							if(x == 1 && y == 0 && Px != 1) // exchange data with right process
							{								
								for(int i = 0; i < local_sizeY; i++)
								{
									U_right_send[i] = U[(y+i)*local_sizeX+(local_sizeX-1)];
									V_right_send[i] = V[(y+i)*local_sizeX+(local_sizeX-1)];
									
								}
								MPI_Send(U_right_send, local_sizeY, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
								MPI_Send(V_right_send, local_sizeY, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_right_import, local_sizeY, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_right_import, local_sizeY, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							
							if(x == local_sizeX-1)
							{
									u_right = U_right_import[y];
									v_right = V_right_import[y];
							}
							else{
									u_right = U[x+1+y*local_sizeX];
									v_right = V[x+1+y*local_sizeX];
							}
							
							if(y == 0 && x == 1) // exchange data with below process
							{								
								for(int i = 0; i < tempXsize-1; i++)
								{
									U_down_send[i] = U[y*local_sizeX+x+i];
									V_down_send[i] = V[y*local_sizeX+x+i];
								}
								
								MPI_Send(U_down_send, tempXsize-1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_down_send, tempXsize-1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_down_import, tempXsize-1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_down_import, tempXsize-1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							
							if(y == 0)
							{
									u_down = U_down_import[x-1];
									v_down = V_down_import[x-1];
							}
							else{
									u_down = U[x+(y-1)*local_sizeX];
									v_down = V[x+(y-1)*local_sizeX];
							}
							
							if(y == 0 && x == 1) // exchange data with above process
							{
								for(int i = 0; i < tempXsize-1; i++)
								{
									U_up_send[i] = U[(local_sizeY-1)*local_sizeX+x+i];
									V_up_send[i] = V[(local_sizeY-1)*local_sizeX+x+i];
								}
								
								MPI_Send(U_up_send, tempXsize-1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_up_send, tempXsize-1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_up_import, tempXsize-1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_up_import, tempXsize-1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							
							if(y == local_sizeY -1)
							{
									u_up = U_up_import[x-1];
									v_up = V_up_import[x-1];
							}
							else{
									u_up = U[x+(y+1)*local_sizeX];
									v_up = V[x+(y+1)*local_sizeX];
							}
							
							U_next[x+y*local_sizeX] = Dt * (
										c*( u_right - 2*U[x+y*local_sizeX] + U[x-1+y*local_sizeX] ) / (Dx*Dx) + 
										c*( u_up - 2*U[x+y*local_sizeX] + u_down ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - U[x-1+y*local_sizeX]) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( v_right - 2*V[x+y*local_sizeX] + V[x-1+y*local_sizeX] ) / (Dx*Dx) + 
										c*( v_up - 2*V[x+y*local_sizeX] + v_down) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - V[x-1+y*local_sizeX]) / Dx 
											  ) + V[x+y*local_sizeX];	
						}
					}
					
					for (int i = 0; i < local_sizeX*local_sizeY; i++)
					{
						U[i] = U_next[i];						
						V[i] = V_next[i];
					}
				}
				
			} // Close processes on the left edge w/o corners
			
		} //Close processes on the left edge
		else if((rank+1)%Px == 0)// Now consider the domain's right edge processes 
		{
			if(rank == Px-1) // right lower corner process
			{
				if(Py == 1)
				{
					tempYsize = local_sizeY-1;
				}
				else
				{
					tempYsize = local_sizeY;
					U_up_send = new double[local_sizeX -1];
					V_up_send = new double[local_sizeX -1];
					U_up_import = new double[local_sizeX -1];
					V_up_import = new double[local_sizeX -1];
				}
				
				U_left_send = new double[tempYsize -1];
				V_left_send = new double[tempYsize -1];
				U_left_import = new double[tempYsize -1];
				V_left_import = new double[tempYsize -1];
				
				for(int T = 1; T < Nt; T++)
				{
					for(int y = 1; y < tempYsize; y++)
					{
						for (int x = 0; x < local_sizeX - 1; x++)
						{
							
							if(x == 0 && y == 1) // exchange data with left process
							{								
								for(int i = 0; i < tempYsize-1; i++)
								{
									U_left_send[i] = U[(y+i)*local_sizeX+x];
									V_left_send[i] = V[(y+i)*local_sizeX+x];
								}
								
								MPI_Send(U_left_send, tempYsize -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
								MPI_Send(V_left_send, tempYsize -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_left_import, tempYsize -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_left_import, tempYsize -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							
							if(x == 0)
							{
									u_left = U_left_import[y-1];
									v_left = V_left_import[y-1];
							}
							else{
									u_left = U[x-1+y*local_sizeX];
									v_left = V[x-1+y*local_sizeX];
							}
							
							if(y == 1 && x == 0 && Py != 1) // exchange data with above process
							{
								for(int i = 0; i < local_sizeX -1; i++)
								{
									U_up_send[i] = U[(local_sizeY-1)*local_sizeX+x+i];
									V_up_send[i] = V[(local_sizeY-1)*local_sizeX+x+i];
									
								}
								
								MPI_Send(U_up_send, local_sizeX -1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_up_send, local_sizeX -1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_up_import, local_sizeX -1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_up_import, local_sizeX -1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							
							}
							
							if(y == local_sizeY -1)
							{
									u_up = U_up_import[x];
									v_up = V_up_import[x];
							}
							else{
									u_up = U[x+(1+y)*local_sizeX];
									v_up = V[x+(1+y)*local_sizeX];
							}
							
							U_next[x+y*local_sizeX] = Dt * (
										c*( U[x+1+y*local_sizeX] - 2*U[x+y*local_sizeX] + u_left ) / (Dx*Dx) + 
										c*( u_up - 2*U[x+y*local_sizeX] + U[x+(y-1)*local_sizeX] ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - U[x+(y-1)*local_sizeX]) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_left) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( V[x+1+y*local_sizeX] - 2*V[x+y*local_sizeX] + v_left ) / (Dx*Dx) + 
										c*( v_up - 2*V[x+y*local_sizeX] + V[x+(y-1)*local_sizeX]) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] - V[x+(y-1)*local_sizeX]) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_left) / Dx 
											  ) + V[x+y*local_sizeX];	
						}
					}
					
					for (int i = 0; i < local_sizeX*local_sizeY; i++)
					{
						U[i] = U_next[i];						
						V[i] = V_next[i];
					}
					
				} 	
				
			} // Close bottom right corner
			else if(rank == (Px*Py-1)) // right top corner
			{
				U_left_send = new double[local_sizeY -1];
				V_left_send = new double[local_sizeY -1];
				U_left_import = new double[local_sizeY -1];
				V_left_import = new double[local_sizeY -1];
				U_down_send = new double[local_sizeX -1];
				V_down_send = new double[local_sizeX -1];
				U_down_import = new double[local_sizeX -1];
				V_down_import = new double[local_sizeX -1];
				
				for(int T = 1; T < Nt; T++)
				{
					for(int y = 0; y < local_sizeY-1; y++)
					{
						for (int x = 0; x < local_sizeX -1 ; x++)
						{
							if(y == 0 && x == 0) // call left adjacent process
							{								
								for(int i = 0; i < local_sizeY -1; i++)
								{
									U_left_send[i] = U[(y+i)*local_sizeX+x];
									V_left_send[i] = V[(y+i)*local_sizeX+x];
								}
								
								MPI_Send(U_left_send, local_sizeY -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
								MPI_Send(V_left_send, local_sizeY -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_left_import, local_sizeY -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_left_import, local_sizeY -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							
							if(x == 0)
							{
									u_left = U_left_import[y];
									v_left = V_left_import[y];
							}
							else{
									u_left = U[x-1+y*local_sizeX];
									v_left = V[x-1+y*local_sizeX];
							}
							
							if(y == 0 && x == 0) // calling down adjacent process
							{								
								for(int i = 0; i < local_sizeX -1; i++)
								{
									U_down_send[i] = U[y*local_sizeX+x+i];
									V_down_send[i] = V[y*local_sizeX+x+i];
									
								}
								
								MPI_Send(U_down_send, local_sizeX -1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_down_send, local_sizeX -1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_down_import, local_sizeX -1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_down_import, local_sizeX -1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							
							}
							
							if(y == 0)
							{
									u_down = U_down_import[x];
									v_down = V_down_import[x];
							}
							else{
									u_down = U[x+(y-1)*local_sizeX];
									v_down = V[x+(y-1)*local_sizeX];
							}
							
							U_next[x+y*local_sizeX] = Dt * (
										c*( U[x+1+y*local_sizeX] - 2*U[x+y*local_sizeX] + u_left ) / (Dx*Dx) + 
										c*( U[x+(y+1)*local_sizeX] - 2*U[x+y*local_sizeX] + u_down ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_left) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( V[x+1+y*local_sizeX] - 2*V[x+y*local_sizeX] + v_left ) / (Dx*Dx) + 
										c*( V[x+(y+1)*local_sizeX] - 2*V[x+y*local_sizeX] + v_down) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_left) / Dx 
											  ) + V[x+y*local_sizeX];	
						}
					}
					
					for (int i = 0; i < local_sizeX*local_sizeY; i++)
					{
						U[i] = U_next[i];						
						V[i] = V_next[i];
					}
				}
				
			} // Close right top corner
			else // Processes on the right side w/o corners
			{
				for(int T = 1; T < Nt; T++)
				{
					for (int y = 0; y < local_sizeY ; y++)
					{
						for(int x = 0; x < local_sizeX-1; x++)
						{
							if(x == 0 && y == 0) // exchange data with left process
							{
								U_left_send = new double[local_sizeY];
								V_left_send = new double[local_sizeY];
								U_left_import = new double[local_sizeY];
								V_left_import = new double[local_sizeY];
								
								for(int i = 0; i < local_sizeY; i++)
								{
									U_left_send[i] = U[(y+i)*local_sizeX+x];
									V_left_send[i] = V[(y+i)*local_sizeX+x];
								}
								
								MPI_Send(U_left_send, local_sizeY, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
								MPI_Send(V_left_send, local_sizeY, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_left_import, local_sizeY, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_left_import, local_sizeY, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							
							}
							if(x == 0)
							{
									u_left = U_left_import[y];
									v_left = V_left_import[y];
							}
							else{
									u_left = U[x-1+y*local_sizeX];
									v_left = V[x-1+y*local_sizeX];
							}
							
							if(y == 0 && x == 0) // calling down adjacent process
							{
								U_down_send = new double[local_sizeX -1];
								V_down_send = new double[local_sizeX -1];
								U_down_import = new double[local_sizeX -1];
								V_down_import = new double[local_sizeX -1];
								
								for(int i = 0; i < local_sizeX -1; i++)
								{
									U_down_send[i] = U[y*local_sizeX+x+i];
									V_down_send[i] = V[y*local_sizeX+x+i];
									
								}
								
								MPI_Send(U_down_send, local_sizeX -1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_down_send, local_sizeX -1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_down_import, local_sizeX -1, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_down_import, local_sizeX -1, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							
							}
							
							if(y == 0)
							{
									u_down = U_down_import[x];
									v_down = V_down_import[x];
							}
							else{
									u_down = U[x+(y-1)*local_sizeX];
									v_down = V[x+(y-1)*local_sizeX];
							}
							
							if(y == 0 && x == 0) // exchange data with above process
							{
								U_up_send = new double[local_sizeX -1];
								V_up_send = new double[local_sizeX -1];
								U_up_import = new double[local_sizeX -1];
								V_up_import = new double[local_sizeX -1];
								
								for(int i = 0; i < local_sizeX -1; i++)
								{
									U_up_send[i] = U[(local_sizeY-1)*local_sizeX+x+i];
									V_up_send[i] = V[(local_sizeY-1)*local_sizeX+x+i];
									
								}
								
								MPI_Send(U_up_send, local_sizeX -1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD);
								MPI_Send(V_up_send, local_sizeX -1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD);
								// MPI_Recv( var, count, type, src,  tag, comm, status)
								MPI_Recv(U_up_import, local_sizeX -1, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Recv(V_up_import, local_sizeX -1, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							
							}
							
							if(y == local_sizeY -1)
							{
									u_up = U_up_import[x];
									v_up = V_up_import[x];
							}
							else{
									u_up = U[x+(1+y)*local_sizeX];
									v_up = V[x+(1+y)*local_sizeX];
							}
							
							U_next[x+y*local_sizeX] = Dt * (
										c*( U[x+1+y*local_sizeX] - 2*U[x+y*local_sizeX] + u_left ) / (Dx*Dx) + 
										c*( u_up - 2*U[x+y*local_sizeX] + u_down ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_left) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( V[x+1+y*local_sizeX] - 2*V[x+y*local_sizeX] + v_left) / (Dx*Dx) + 
										c*( v_up - 2*V[x+y*local_sizeX] + v_down) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] -v_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_left) / Dx 
											  ) + V[x+y*local_sizeX];
							
						}
					}
					
					for (int i = 0; i < local_sizeX*local_sizeY; i++)
					{
						U[i] = U_next[i];						
						V[i] = V_next[i];
					}
					
				}
				
			} // Close right edge processes w/o corners
			
		} // Close processes on the right edge of domain  
		else if(rank > 0 && rank < Px-1) // Now consider bottom edge processes, w/o corners
		{
			if(Py == 1)
			{
				tempYsize = local_sizeY-1;
			}
			else
			{
				tempYsize = local_sizeY;
				U_up_send = new double[local_sizeX];
				V_up_send = new double[local_sizeX];
				U_up_import = new double[local_sizeX];
				V_up_import = new double[local_sizeX];
			}
			
			U_left_send = new double[tempYsize -1];
			V_left_send = new double[tempYsize -1];
			U_left_import = new double[tempYsize -1];
			V_left_import = new double[tempYsize -1];
			U_right_send = new double[tempYsize -1];
			V_right_send = new double[tempYsize -1];
			U_right_import = new double[tempYsize -1];
			V_right_import = new double[tempYsize -1];				
							
			for(int T = 1; T < Nt; T++)
			{
				for (int y = 1; y < tempYsize ; y++)
				{
					for(int x = 0; x < local_sizeX; x++)
					{
						
						if(y == 1 && x == 0 && Py != 1) // calling above adjacent process
						{
							for(int i = 0; i < local_sizeX ; i++)
							{
								U_up_send[i] = U[(local_sizeY-1)*local_sizeX+x+i];
								V_up_send[i] = V[(local_sizeY-1)*local_sizeX+x+i];
							}
							
							MPI_Send(U_up_send, local_sizeX, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD);
							MPI_Send(V_up_send, local_sizeX, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_up_import, local_sizeX, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_up_import, local_sizeX, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						
						if(y == local_sizeY-1)
						{
								u_up = U_up_import[x];
								v_up = V_up_import[x];
						}
						else{
								u_up = U[x+(y+1)*local_sizeX];
								v_up = V[x+(y+1)*local_sizeX];
						}
						
						if(y == 1 && x == 0) // calling left adjacent process
						{
							for(int i = 0; i < tempYsize -1; i++)
							{
								U_left_send[i] = U[(y+i)*local_sizeX+x];
								V_left_send[i] = V[(y+i)*local_sizeX+x];
							}
							
							MPI_Send(U_left_send, tempYsize -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
							MPI_Send(V_left_send, tempYsize -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_left_import, tempYsize -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_left_import, tempYsize -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						
						if(x == 0)
						{
								u_left = U_left_import[y-1];
								v_left = V_left_import[y-1];
						}
						else{
								u_left = U[x-1+y*local_sizeX];
								v_left = V[x-1+y*local_sizeX];
						}
						
						if(y == 1 && x == 0) // calling right adjacent process
						{
							for(int i = 0; i < tempYsize -1; i++)
							{
								U_right_send[i] = U[(y+i)*local_sizeX+local_sizeX-1];
								V_right_send[i] = V[(y+i)*local_sizeX+local_sizeX-1];
							}
							
							MPI_Send(U_right_send, tempYsize -1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
							MPI_Send(V_right_send, tempYsize -1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_right_import, tempYsize -1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_right_import, tempYsize -1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						
						}
						
						if(x == local_sizeX-1)
						{
								u_right = U_right_import[y-1];
								v_right = V_right_import[y-1];
						}
						else{
								u_right = U[x+1+y*local_sizeX];
								v_right = V[x+1+y*local_sizeX];
						}
						
							U_next[x+y*local_sizeX] = Dt * (
										c*( u_right - 2*U[x+y*local_sizeX] + u_left ) / (Dx*Dx) + 
										c*( u_up - 2*U[x+y*local_sizeX] + U[x+(y-1)*local_sizeX] ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - U[x+(y-1)*local_sizeX]) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_left) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( v_right - 2*V[x+y*local_sizeX] + v_left) / (Dx*Dx) + 
										c*( v_up - 2*V[x+y*local_sizeX] + V[x+(y-1)*local_sizeX]) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] -V[x+(y-1)*local_sizeX]) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_left) / Dx 
											  ) + V[x+y*local_sizeX];
						
					}
				}
					
				for (int i = 0; i < local_sizeX*local_sizeY; i++)
				{
					U[i] = U_next[i];						
					V[i] = V_next[i];
				}
				
			}
			
		} // Close lower edge processes w/o corners
		else if(rank > (Py-1)*Px && rank < Py*Px-1 ) // top edge processes w/o corners processes
		{
			U_down_send = new double[local_sizeX];
			V_down_send = new double[local_sizeX];
			U_down_import = new double[local_sizeX];
			V_down_import = new double[local_sizeX];
			U_left_send = new double[local_sizeY -1];
			V_left_send = new double[local_sizeY -1];
			U_left_import = new double[local_sizeY -1];
			V_left_import = new double[local_sizeY -1];
			U_right_send = new double[local_sizeY -1];
			V_right_send = new double[local_sizeY -1];
			U_right_import = new double[local_sizeY -1];
			V_right_import = new double[local_sizeY -1];
							
			for(int T = 1; T < Nt; T++)
			{
				for (int y = 0; y < local_sizeY-1 ; y++)
				{
					for(int x = 0; x < local_sizeX; x++)
					{
						if(y == 0 && x == 0) // calling below adjacent process
						{
							for(int i = 0; i < local_sizeX ; i++)
							{
								U_down_send[i] = U[y*local_sizeX+x+i];
								V_down_send[i] = V[y*local_sizeX+x+i];
							}
							
							MPI_Send(U_down_send, local_sizeX, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD);
							MPI_Send(V_down_send, local_sizeX, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_down_import, local_sizeX, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_down_import, local_sizeX, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						
						if(y == 0)
						{
								u_down = U_down_import[x];
								v_down = V_down_import[x];
						}
						else{
								u_down = U[x+(y-1)*local_sizeX];
								v_down = V[x+(y-1)*local_sizeX];
						}
						
						if(y == 0 && x == 0) // calling left adjacent process
						{
							for(int i = 0; i < local_sizeY -1; i++)
							{
								U_left_send[i] = U[(y+i)*local_sizeX+x];
								V_left_send[i] = V[(y+i)*local_sizeX+x];
							}
							
							MPI_Send(U_left_send, local_sizeY -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
							MPI_Send(V_left_send, local_sizeY -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_left_import, local_sizeY -1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_left_import, local_sizeY -1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						
						}
						
						if(x == 0)
						{
								u_left = U_left_import[y];
								v_left = V_left_import[y];
						}
						else{
								u_left = U[x-1+y*local_sizeX];
								v_left = V[x-1+y*local_sizeX];
						}
						
						
						if(y == 0 && x == 0) // calling right adjacent process
						{							
							for(int i = 0; i < local_sizeY -1; i++)
							{
								U_right_send[i] = U[(y+i)*local_sizeX+local_sizeX-1];
								V_right_send[i] = V[(y+i)*local_sizeX+local_sizeX-1];
								
							}
							
							MPI_Send(U_right_send, local_sizeY -1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
							MPI_Send(V_right_send, local_sizeY -1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_right_import, local_sizeY -1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_right_import, local_sizeY -1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						
						}
						
						if(x == local_sizeX-1)
						{
								u_right = U_right_import[y];
								v_right = V_right_import[y];
						}
						else{
								u_right = U[x+1+y*local_sizeX];
								v_right = V[x+1+y*local_sizeX];
						}
						
							U_next[x+y*local_sizeX] = Dt * (
										c*( u_right - 2*U[x+y*local_sizeX] + u_left ) / (Dx*Dx) + 
										c*( U[x+(y+1)*local_sizeX] - 2*U[x+y*local_sizeX] + u_down ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_left) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( v_right - 2*V[x+y*local_sizeX] + v_left) / (Dx*Dx) + 
										c*( V[x+(y+1)*local_sizeX] - 2*V[x+y*local_sizeX] + v_down) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] -v_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_left) / Dx 
											  ) + V[x+y*local_sizeX];
						
						// Integrate
					}
				}
					
				for (int i = 0; i < local_sizeX*local_sizeY; i++)
				{
					U[i] = U_next[i];						
					V[i] = V_next[i];
				}
				
			}

		} // Close ranks whoch belong to the top edge w/o corners
		else // At last consider the middle processes, those inside domain (not boundary) 
		{
			U_right_send = new double[local_sizeY];
			V_right_send = new double[local_sizeY];
			U_right_import = new double[local_sizeY];
			V_right_import = new double[local_sizeY];
			U_left_send = new double[local_sizeY];
			V_left_send = new double[local_sizeY];
			U_left_import = new double[local_sizeY];
			V_left_import = new double[local_sizeY];
			U_down_send = new double[local_sizeX];
			V_down_send = new double[local_sizeX];
			U_down_import = new double[local_sizeX];
			V_down_import = new double[local_sizeX];
			U_up_send = new double[local_sizeX];
			V_up_send = new double[local_sizeX];
			U_up_import = new double[local_sizeX];
			V_up_import = new double[local_sizeX];
								
			for(int T = 1; T < Nt; T++)
			{
				for (int y = 0; y < local_sizeY ; y++)
				{
					for(int x = 0; x < local_sizeX; x++)
					{
						
						if(y == 0 && x == 0) // calling right adjacent process
						{
							for(int i = 0; i < local_sizeY; i++)
							{
								U_right_send[i] = U[(y+i)*local_sizeX+local_sizeX-1];
								V_right_send[i] = V[(y+i)*local_sizeX+local_sizeX-1];
							}
							
							MPI_Send(U_right_send, local_sizeY , MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
							MPI_Send(V_right_send, local_sizeY , MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_right_import, local_sizeY, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_right_import, local_sizeY, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						
						if(x == local_sizeX-1)
						{
								u_right = U_right_import[y];
								v_right = V_right_import[y];
						}
						else{
								u_right = U[x+1+y*local_sizeX];
								v_right = V[x+1+y*local_sizeX];
						} 
						
						
						if(y == 0 && x == 0) // calling left adjacent process
						{							
							for(int i = 0; i < local_sizeY; i++)
							{
								U_left_send[i] = U[(y+i)*local_sizeX+x];
								V_left_send[i] = V[(y+i)*local_sizeX+x];
							}
							
							MPI_Send(U_left_send, local_sizeY, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
							MPI_Send(V_left_send, local_sizeY, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_left_import, local_sizeY, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_left_import, local_sizeY, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						
						}
						
						if(x == 0)
						{
								u_left = U_left_import[y];
								v_left = V_left_import[y];
						}
						else{
								u_left = U[x-1+y*local_sizeX];
								v_left = V[x-1+y*local_sizeX];
						} 
						
						if(y == 0 && x == 0) // calling below adjacent process
						{							
							for(int i = 0; i < local_sizeX ; i++)
							{
								U_down_send[i] = U[y*local_sizeX+x+i];
								V_down_send[i] = V[y*local_sizeX+x+i];
							}
							
							MPI_Send(U_down_send, local_sizeX, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD);
							MPI_Send(V_down_send, local_sizeX, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_down_import, local_sizeX, MPI_DOUBLE, rank-Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_down_import, local_sizeX, MPI_DOUBLE, rank-Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						
						if(y == 0)
						{
								u_down = U_down_import[x];
								v_down = V_down_import[x];
						}
						else{
								u_down = U[x+(y-1)*local_sizeX];
								v_down = V[x+(y-1)*local_sizeX];
						}
						
						if(y == 0 && x == 0) // exchange data with above process
						{
							for(int i = 0; i < local_sizeX; i++)
							{
								U_up_send[i] = U[(local_sizeY-1)*local_sizeX+x+i];
								V_up_send[i] = V[(local_sizeY-1)*local_sizeX+x+i];
							}
							
							MPI_Send(U_up_send, local_sizeX, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD);
							MPI_Send(V_up_send, local_sizeX, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD);
							// MPI_Recv( var, count, type, src,  tag, comm, status)
							MPI_Recv(U_up_import, local_sizeX, MPI_DOUBLE, rank+Px, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Recv(V_up_import, local_sizeX, MPI_DOUBLE, rank+Px, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
							
							if(y == local_sizeY -1)
							{
									u_up = U_up_import[x];
									v_up = V_up_import[x];
							}
							else{
									u_up = U[x+(1+y)*local_sizeX];
									v_up = V[x+(1+y)*local_sizeX];
							}
							
							U_next[x+y*local_sizeX] = Dt * (
										c*( u_right - 2*U[x+y*local_sizeX] + u_left ) / (Dx*Dx) + 
										c*( u_up - 2*U[x+y*local_sizeX] + u_down ) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (U[x+y*local_sizeX] - u_left) / Dx 
											  ) + U[x+y*local_sizeX];
										
							V_next[x+y*local_sizeX] = Dt * ( 
										c*( v_right - 2*V[x+y*local_sizeX] + v_left) / (Dx*Dx) + 
										c*( v_up - 2*V[x+y*local_sizeX] + v_down) / (Dy*Dy) -
										(ay + b*V[x+y*local_sizeX]) * (V[x+y*local_sizeX] -v_down) / Dy -
										(ax + b*U[x+y*local_sizeX]) * (V[x+y*local_sizeX] - v_left) / Dx 
											  ) + V[x+y*local_sizeX];
					}
				}
					
				for (int i = 0; i < local_sizeX*local_sizeY; i++)
				{
					U[i] = U_next[i];						
					V[i] = V_next[i];
				}
				
			} // Close integration Loop
		
		} // Close domain rank scanning
		delete[] U_next;
		delete[] V_next;
	
	} // Close multi-process integration
		
} // Close integration function 

void Burgers::WriteToFile(int rank, int size)
{
	if(size == 1) //  This method is tailored for single-process, as it omits unnecessary MPI_Send and MPI_Recv
	{
	
		ofstream plik("VelocityField.txt", ios::out| ios::trunc);
		if (!plik.is_open())
		{
			cout << "Failed to open an output file" << endl;
		}
		else
		{	
			if (plik.good())
			{
				plik << setw(16) << "X" << setw(16) << "Y" << setw(16) << "U" << setw(16) << "V" << endl;
				plik << endl; // make space in the output to improve readibility
			
				for (int j = 0; j < sizeY; j++)
				{	
					for(int i = 0; i < sizeX; i++)
					{
						plik.precision(4);
						plik << setw(16) << X_local[i] << setw(16) << Y_local[j] << setw(16) << U[i+j*sizeX] << setw(16) << V[i+j*sizeX] << endl;
					
					}
					plik << endl; 
				}	
				plik.close();
			}
			
		}
		
	}
	else// For multiple-process computation
	{
		// Each process except 0, sends its arrays to the root process
		if(rank != 0)
		{
		MPI_Send(&local_sizeX, 1 , MPI_INT, 0 ,0, MPI_COMM_WORLD);
		MPI_Send(&local_sizeY, 1 , MPI_INT, 0 ,1, MPI_COMM_WORLD);
		MPI_Send(U,local_sizeX*local_sizeY , MPI_DOUBLE, 0 ,2, MPI_COMM_WORLD);
		MPI_Send(V,local_sizeX*local_sizeY , MPI_DOUBLE, 0 ,3, MPI_COMM_WORLD);
		MPI_Send(X_local, local_sizeX , MPI_DOUBLE, 0 ,4, MPI_COMM_WORLD);
		MPI_Send(Y_local, local_sizeY , MPI_DOUBLE, 0 ,5, MPI_COMM_WORLD);
		}
			
		if(rank == 0)
		{
			ofstream plik("VelocityField.txt", ios::out| ios::trunc);
			if (!plik.is_open())
			{
				cout << "Failed to open an output file" << endl;
			}
			else
			{	
				if (plik.good())
				{
					//plik << setw(16) << "X" << setw(16) << "Y" << setw(16) << "U" << setw(16) << "V" << endl;
					plik << endl;
					
					for (int j = 0; j < local_sizeY; j++)
					{	
						for(int i = 0; i < local_sizeX; i++)
						{
							plik.precision(4);
							plik << setw(16) << fixed << X_local[i] << setw(16) << fixed<< Y_local[j] << setw(16)<< fixed 
								 << U[i+j*local_sizeX] << setw(16) << fixed << V[i+j*local_sizeX] << endl;
						}
						
						plik << endl;
					}
					
				}
				
				for(int p = 1; p < size; p++)
				{	
					MPI_Recv(&import_xSize, 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);				
					MPI_Recv(&import_ySize, 1, MPI_INT, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					U_import = new double[import_xSize*import_ySize];
					V_import = new double[import_xSize*import_ySize];
					X_import = new double[import_xSize];
					Y_import = new double[import_ySize];
					MPI_Recv(U_import, import_xSize*import_ySize , MPI_DOUBLE, p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(V_import, import_xSize*import_ySize , MPI_DOUBLE, p, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(X_import, import_xSize , MPI_DOUBLE, p, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(Y_import, import_ySize , MPI_DOUBLE, p, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					
					for (int j = 0; j < import_ySize; j++)
					{	
						for(int i = 0; i < import_xSize; i++)
						{
							if (plik.good())
							{
								plik.precision(4);
								plik << setw(16) << fixed << X_import[i] << setw(16) << fixed << Y_import[j] << setw(16)<< fixed 
									 << U_import[i+j*import_xSize] << setw(16) << fixed<< V_import[i+j*import_xSize] << endl;
							}
							
						}
						plik << endl; // make space in the output to improve readibility
					}
					// Free memory and prepare for next iteeration
					delete[] X_import;
					delete[] Y_import;
					delete[] U_import;
					delete[] V_import;
				}
				
				plik.close();
			}
		}	
	}
}

void Burgers::CalcEnergy(int rank, int size)
{
	if(size == 1) // This method is for single-process as it omits MPI, hence is quicker
	{
		// This is the quickest way to get the energy, explained in the report
		ALLuNorm = cblas_dnrm2(sizeXY,U,1); 
		ALLvNorm = cblas_dnrm2(sizeXY,V,1);
		Total_Energy = 0.5 * (pow(ALLuNorm,2) + pow(ALLvNorm,2)) *Dx*Dy;
		cout << "Total energy is: " << Total_Energy << endl;
	}
	else // This method is for multi-process computation
	{
		ALLuNorm = cblas_dnrm2(local_sizeX*local_sizeY,U,1); 
		ALLvNorm = cblas_dnrm2(local_sizeX*local_sizeY,V,1);
		Energy = 0.5 * (pow(ALLuNorm,2) + pow(ALLvNorm,2)) *Dx*Dy;
		MPI_Reduce(&Energy, &Total_Energy, 1, MPI_DOUBLE, MPI_SUM , 0, MPI_COMM_WORLD);
			
		if (rank == 0)
		{
			cout << endl << "Total energy is: " << Total_Energy << endl;
		}
	}
	
}
