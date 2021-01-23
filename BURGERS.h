#ifndef CLASS_BURGERS
#define CLASS_BURGERS
#include "MODEL.h"

class Burgers
{
private:

	double *U;
	double *V;
	double *U_next;
	double *V_next;
	
	double ALLuNorm;
	double ALLvNorm;
	
	double *U_import;
	double *V_import;
	double *Y_import;
	double *X_import;
	int import_xSize;
	int import_ySize;
	
	double Dx;
	double Dy;
	double Dt;
	
	int Nt;
	int Px;
	int Py;
	int sizeX;
	int sizeY;
	int sizeXY;
	int local_sizeX;
	int local_sizeY;
	double* Y_local;
	double* X_local;
	int tempXsize;
	int tempYsize;
	
	double* U_right_import;
	double* V_right_import;
	double* U_left_import;
	double* V_left_import;
	double* U_up_import;
	double* V_up_import;
	double* U_down_import;
	double* V_down_import;
	
	double* U_right_send;
	double* V_right_send;
	double* U_left_send;
	double* V_left_send;
	double* U_up_send;
	double* V_up_send;
	double* U_down_send;
	double* V_down_send;
	
	double u_right;
	double v_right;
	double u_left;
	double v_left;
	double u_up;
	double v_up;
	double u_down;
	double v_down;
	
	
	double Lx_half;
	double Ly_half;
	
	double max_xID;
	double min_xID;
	double max_yID;
	double min_yID;
	
	double radius;
	double ax;
	double ay;
	double b;
	double c;
	double Energy;
	double Total_Energy;
	
public:
	Burgers(const Model& model);
	Burgers();
   ~Burgers();
	void VelocityVectors(int);
	void SpaceVectors(int, int);
	void Integrate(int ,int);
	void WriteToFile(int, int);
	void CalcEnergy(int, int);	
};

#endif