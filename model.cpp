#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <ctime>
#include <stdio.h>
#include <time.h>
#include <filesystem>
using namespace std;

clock_t t_0;
clock_t t_1;

bool check_time_waste = false;
bool save_in_xyz = true;

double delta_t = 1/1000.0;

double r_consid = 5; //epsilon approx 0.15%

const string currentDateTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

class atom 
{
private:
	bool Real; //is atom from real box or from non-existing box
	double r[3]; //cord
	double v[3]; //velocity
	double a[3]; //acceleration
	double E; //Kinetic Energy
	double delta[3]; //Vect between atoms
	double d; //distance between atoms
	double d_squared;//d^2
	double U; //Potential Energy
	double F_val; //F between atoms
	double box_size; //working box
public:
atom() {}
atom(double (&r0)[3], double (&v0)[3], double (&a0)[3], double box_size0, bool Real0)
	{
		r[0] = r0[0];
		r[1] = r0[1];
		r[2] = r0[2];
		v[0] = v0[0];
		v[1] = v0[1];
		v[2] = v0[2];
		a[0] = a0[0];
		a[1] = a0[1];
		a[2] = a0[2];
		box_size = box_size0;
		Real = Real0;
	E = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0;
	}

void Make_Zero_a()
{
	r[0] = 0;
	r[1] = 0;
	r[2] = 0;
}

double* Get_r() { return r; }
double* Get_v() { return v; }
double* Get_a() { return a; }
double Get_E() { return E; }
double Get_U() {return U; }

void Set_r(double (&r0)[3]){ 
r[0] = r0[0];
r[1] = r0[1];
r[2] = r0[2];
if (Real)
{
if (r[0] <= 0)
{
r[0] = box_size + r[0];
}
if (r[0] > box_size)
{
r[0] = r[0] - box_size;
}
if (r[1] <= 0)
{
r[1] = box_size + r[1];
}
if (r[1] > box_size)
{
r[1] = r[1] - box_size;
}
if (r[2] <= 0)
{
r[2] = box_size + r[2];
}
if (r[2] > box_size)
{
r[2] = r[2] - box_size;
}
}
}
void Set_v(double (&v0)[3]) {
v[0] = v0[0];
v[1] = v0[1];
v[2] = v0[2];
}
void Set_a(double (&a0)[3]) {
a[0] = a0[0];
a[1] = a0[1];
a[2] = a0[2];
}

void Set_U(double U0) {
U = U0;
}

void Update_E()
{
	E = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0;
}


void Calculate_F(atom temp_atom, double d_squared)
{
	delta[0] = r[0] - temp_atom.r[0];
	delta[1] = r[1] - temp_atom.r[1];
	delta[2] = r[2] - temp_atom.r[2];
	if (d_squared==0) {d=0;} else {d=pow(d_squared, 0.5);}
	F_val = -4*(6/(d*d*d*d*d*d*d)-12/(d*d*d*d*d*d*d*d*d*d*d*d*d));
	//F = -4 epsilon ( 6 a^6/r^7 - 12 a^12/r^13 )
	a[0] += F_val*delta[0]/d;
	a[1] += F_val*delta[1]/d;
	a[2] += F_val*delta[2]/d;
	//cout << delta[0] << " " << delta[1] << " " << delta[2] << endl;
}
void Calculate_U(double d_squared) {
	if (d_squared==0) {d=0;} else {d=pow(d_squared, 0.5);}
	U += 4*(-1/(d*d*d*d*d*d)+1/(d*d*d*d*d*d*d*d*d*d*d*d));
}
void All_Info()
{
	cout << "r: " << r[0] << " " << r[1] << " " << r[2] << endl;
	cout << "v: " << v[0] << " " << v[1] << " " << v[2] << endl;
	cout << "a: " << a[0] << " " << a[1] << " " << a[2] << endl;
}
};
class Group_Of_Atoms
{
	private: 
		//погуглить про связанные листы
		vector<atom> atoms;
		vector<atom> fantom_atoms;
		string dir_name;
		double E;
		double U;
		double temp_r[3];
		double temp_v[3];
		double temp_a[3];
		int length;
		double d_temp;
		double box_size;
	public:
		Group_Of_Atoms(int n, string filename)
		{
		dir_name = filename;
		box_size = 3*(n-1)+1;
		double pos[3]{0, 0, 0};
		double zero_vec[3]{0, 0, 0};
			for (int i=0; i<n; i++)
			{
				for (int j=0; j<n; j++)
				{
					for (int k=0; k<n; k++)
					{
						pos[0] = 3*i + 0.5;
						pos[1] = 3*j + 0.5;
						pos[2] = 3*k + 0.5;
					atom temporary_atom(pos, zero_vec, zero_vec, box_size, true);
					atoms.push_back(temporary_atom);
					}
				}
			}
		length = atoms.size();
		}

	
		void Group_Info()
		{
			for (int i=0; i<length; i++)
			{
				atoms[i].All_Info();
			}
		}
		atom Get_Atom(int index)
		{
			return atoms[index];
		}
		void Fantom_Group_Info()
		{
			for (int i=0; i<fantom_atoms.size(); i++)
			{
				fantom_atoms[i].All_Info();
			}
		}
		void Update_Fantom_Atoms()
		{
			fantom_atoms.clear();
			for (int i=0; i<length; i++)
			{
				temp_r[0] = atoms[i].Get_r()[0];
				temp_r[1] = atoms[i].Get_r()[1];
				temp_r[2] = atoms[i].Get_r()[2];
				for (double a : {-1, 0, 1})
				{
					for (double b : {-1, 0, 1})
					{
						for (double c : {-1, 0, 1})
						{
							if (a*a + b*b + c*c != 0)
							{
							temp_r[0] = atoms[i].Get_r()[0] + box_size*a;
							temp_r[1] = atoms[i].Get_r()[1] + box_size*b;
							temp_r[2] = atoms[i].Get_r()[2] + box_size*c;
							//cout << a << " " << b << " " << c << endl;
							
							if ((temp_r[0] <= box_size+r_consid)&&(temp_r[0] > -1*r_consid)&&(temp_r[1] <= box_size+r_consid)&&(temp_r[1] > -1*r_consid)&&(temp_r[2] <= box_size+r_consid)&&(temp_r[2] > -1*r_consid))
							{
								atom temporary_atom(temp_r, temp_r, temp_r, box_size, false);
								fantom_atoms.push_back(temporary_atom);
							}
							}
						}
					}
				}
			}
			//Fantom_Group_Info();
		}
		void One_Velocity_Verlet_Iteration()
		{
			E = 0;
			U = 0;
			vector<vector<double>> temp;
			t_0 = clock();
			for (int i=0; i<length; i++) 
			{
				temp_v[0] = atoms[i].Get_v()[0]+delta_t*(1/2.0)*atoms[i].Get_a()[0];
				temp_v[1] = atoms[i].Get_v()[1]+delta_t*(1/2.0)*atoms[i].Get_a()[1];
				temp_v[2] = atoms[i].Get_v()[2]+delta_t*(1/2.0)*atoms[i].Get_a()[2];
				atoms[i].Set_v(temp_v);
			}
			if (check_time_waste)
			{
				cout << "1:	" << (double)(clock() - t_0)/CLOCKS_PER_SEC << endl;
			}
			t_0 = clock();
			for (int i=0; i<length; i++) 
			{
				temp_r[0] = atoms[i].Get_r()[0] + delta_t*(1/2.0)*atoms[i].Get_v()[0];
				temp_r[1] = atoms[i].Get_r()[1] + delta_t*(1/2.0)*atoms[i].Get_v()[1];
				temp_r[2] = atoms[i].Get_r()[2] + delta_t*(1/2.0)*atoms[i].Get_v()[2];
				atoms[i].Set_r(temp_r);
			}
			if (check_time_waste)
			{
				cout << "2:	" << (double)(clock() - t_0)/CLOCKS_PER_SEC << endl;
			}
			t_0 = clock();
			for (int i=0; i<length; i++)
			{
				temp_a[0] = 0;
				temp_a[1] = 0;
				temp_a[2] = 0;
				atoms[i].Set_a(temp_a);
			}
			Update_Fantom_Atoms();
			for (int i=0; i<length; i++) 
			{
				t_1 = clock();
				for (int j = 0; j<length; j++)
				{
					d_temp = (atoms[i].Get_r()[0] - atoms[j].Get_r()[0])*(atoms[i].Get_r()[0] - atoms[j].Get_r()[0]) +  (atoms[i].Get_r()[1] - atoms[j].Get_r()[1])*(atoms[i].Get_r()[1] - atoms[j].Get_r()[1]) + (atoms[i].Get_r()[2] - atoms[j].Get_r()[2])*(atoms[i].Get_r()[2] - atoms[j].Get_r()[2]);	
					if ((d_temp!=0)&&(d_temp<=r_consid*r_consid))
					//if ((d_temp!=0))
					{
						atoms[i].Calculate_F(atoms[j], d_temp);
					}
				}
				for (int j = 0; j<fantom_atoms.size(); j++)
				{
					//cout << fantom_atoms[j].Get_r()[0] << endl;
					//fantom_atoms[j].All_Info();
					d_temp = (atoms[i].Get_r()[0] - fantom_atoms[j].Get_r()[0])*(atoms[i].Get_r()[0] - fantom_atoms[j].Get_r()[0]) +  (atoms[i].Get_r()[1] - fantom_atoms[j].Get_r()[1])*(atoms[i].Get_r()[1] - fantom_atoms[j].Get_r()[1]) + (atoms[i].Get_r()[2] - fantom_atoms[j].Get_r()[2])*(atoms[i].Get_r()[2] - fantom_atoms[j].Get_r()[2]);
					//cout << d_temp << endl;
					if ((d_temp!=0)&&(d_temp<=r_consid*r_consid))
					//if ((d_temp!=0))
					{
						atoms[i].Calculate_F(fantom_atoms[j], d_temp);
					}
				}
				if (check_time_waste)
				{
					cout << "3_" << i << ":	" << (double)(clock() - t_1)/CLOCKS_PER_SEC << endl;
				}
			}
			if (check_time_waste)
			{
				cout << "3_total:	" << (double)(clock() - t_0)/CLOCKS_PER_SEC << endl;
			}
			t_0 = clock();
			for (int i=0; i<length; i++) 
			{
				temp_v[0] = atoms[i].Get_v()[0]+delta_t*(1/2.0)*atoms[i].Get_a()[0];
				temp_v[1] = atoms[i].Get_v()[1]+delta_t*(1/2.0)*atoms[i].Get_a()[1];
				temp_v[2] = atoms[i].Get_v()[2]+delta_t*(1/2.0)*atoms[i].Get_a()[2];
				atoms[i].Set_v(temp_v);
			}

			if (check_time_waste)
			{
				cout << "4:	" << (double)(clock() - t_0)/CLOCKS_PER_SEC << endl;
			}
		}
		void Calculate_Energy()
		{
			E = 0;
			U = 0;
			for (int i=0; i<length; i++) 
			{
				atoms[i].Update_E();
				E += atoms[i].Get_E();
				atoms[i].Set_U(0);
				for (int j = 0; j<length; j++)
				{
					d_temp = (atoms[i].Get_r()[0] - atoms[j].Get_r()[0])*(atoms[i].Get_r()[0] - atoms[j].Get_r()[0]) +  (atoms[i].Get_r()[1] - atoms[j].Get_r()[1])*(atoms[i].Get_r()[1] - atoms[j].Get_r()[1]) + (atoms[i].Get_r()[2] - atoms[j].Get_r()[2])*(atoms[i].Get_r()[2] - atoms[j].Get_r()[2]);	
					if ((d_temp!=0)&&(d_temp<=r_consid*r_consid))
					//if ((d_temp!=0))
					{
						atoms[i].Calculate_U(d_temp);
					}
				}
				for (int j = 0; j<fantom_atoms.size(); j++)
				{
					//cout << fantom_atoms[j].Get_r()[0] << endl;
					d_temp = (atoms[i].Get_r()[0] - fantom_atoms[j].Get_r()[0])*(atoms[i].Get_r()[0] - fantom_atoms[j].Get_r()[0]) +  (atoms[i].Get_r()[1] - fantom_atoms[j].Get_r()[1])*(atoms[i].Get_r()[1] - fantom_atoms[j].Get_r()[1]) + (atoms[i].Get_r()[2] - fantom_atoms[j].Get_r()[2])*(atoms[i].Get_r()[2] - fantom_atoms[j].Get_r()[2]);
					if ((d_temp!=0)&&(d_temp<=r_consid*r_consid))
					//if ((d_temp!=0))
					{
						atoms[i].Calculate_U(d_temp);
					}
				}
				U+=atoms[i].Get_U();
			}

		}

		void Write_r_data()
		{
			std::ofstream stream;                    
			stream.open(dir_name + "//r_data.xyz", std::ios::app);  
			stream << to_string(length) << std::endl;
			stream << std::endl;
			for (int i=0; i<length; i++)
			{
				temp_r[0] = atoms[i].Get_r()[0];
				temp_r[1] = atoms[i].Get_r()[1];
				temp_r[2] = atoms[i].Get_r()[2];
				stream << "H	" << temp_r[0] << "	" << temp_r[1] << "	" << temp_r[2] << "	" << endl;
			}
			stream.close();

		}
		void Write_v_data()
		{
			std::ofstream stream;                    
			stream.open(dir_name + "//v_data.xyz", std::ios::app);  
			stream << to_string(length) << std::endl;
			stream << std::endl;
			for (int i=0; i<length; i++)
			{
				temp_v[0] = atoms[i].Get_v()[0];
				temp_v[1] = atoms[i].Get_v()[1];
				temp_v[2] = atoms[i].Get_v()[2];
				stream << "H	" << temp_v[0] << "	" << temp_v[1] << "	" << temp_v[2] << "	" << endl;
			}
			stream.close();

		}

		void Write_Energy_data(int i)
		{	
			std::ofstream stream;                    
			stream.open(dir_name + "//Energy_data.xyz", std::ios::app);  
			stream << i << "	" << U+E << "	" << E << "	" << U << "	" << endl;		
			stream.close();
		}

		void Calculation_And_Writing_Cycle(int time)
		{
			if (save_in_xyz)
			{
				std::ofstream outfile1 (dir_name + "//r_data.xyz");
				outfile1.close();
				std::ofstream outfile2 (dir_name + "//Energy_data.xyz");
				outfile2 << "i_amount	"<< "Total energy	" << "Kinetic energy	" << "Potential energy	" << endl;
				outfile2.close();
				std::ofstream outfile3 (dir_name + "//v_data.xyz");	
				outfile3.close();
				Write_r_data();
				Write_v_data();
			}
			int m = time/delta_t;
			//Group_Info();
			for (int i=0; i<m; i++)
			{
				One_Velocity_Verlet_Iteration();
				Calculate_Energy();
				if (i%100 ==0)
				{
				cout << i << endl;
				//Group_Info();
				if (save_in_xyz)
				{
					Write_r_data();
					Write_v_data();
					Write_Energy_data(i);
				}
				}
			}
		}

};

int main() 
{
	string date = currentDateTime();
	if (save_in_xyz)
	{
	cout << date << endl;
	std::__fs::filesystem::create_directory(date);
	}
	Group_Of_Atoms Group_1(10, date);
	Group_1.Calculation_And_Writing_Cycle(120);
}
