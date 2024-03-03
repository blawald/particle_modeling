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

float delta_t = 1/10000.0;

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
	float r[3];
	float v[3];
	float a[3];
	float E;
	float delta[3];
	float d;
	float d_squared;
	float U;
	float F_val;
public:
atom() {}
atom(float (&r0)[3], float (&v0)[3], float (&a0)[3])
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
	E = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0;
	}

void Make_Zero_a()
{
	r[0] = 0;
	r[1] = 0;
	r[2] = 0;
}

float* Get_r() { return r; }
float* Get_v() { return v; }
float* Get_a() { return a; }
float Get_E() { return E; }
float Get_U() {return U; }

void Set_r(float (&r0)[3]){ 
r[0] = r0[0];
r[1] = r0[1];
r[2] = r0[2];
}
void Set_v(float (&v0)[3]) {
v[0] = v0[0];
v[1] = v0[1];
v[2] = v0[2];
}
void Set_a(float (&a0)[3]) {
a[0] = a0[0];
a[1] = a0[1];
a[2] = a0[2];
}

void Set_U(float U0) {
U = U0;
}

void Update_E()
{
	E = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0;
}


void Calculate_F(atom temp_atom)
{
	delta[0] = r[0] - temp_atom.r[0];
	delta[1] = r[1] - temp_atom.r[1];
	delta[2] = r[2] - temp_atom.r[2];
	d_squared = delta[1]*delta[1] + delta[2]*delta[2] + delta[0]*delta[0];
	if (d_squared==0) {d=0;} else {d=pow(d_squared, 0.5);}
	F_val = -4*(6/(d*d*d*d*d*d*d)-12/(d*d*d*d*d*d*d*d*d*d*d*d*d));
	//F = -4 epsilon ( 6 a^6/r^7 - 12 a^12/r^13 )
	a[0] += F_val*delta[0]/d;
	a[1] += F_val*delta[1]/d;
	a[2] += F_val*delta[2]/d;
	//cout << delta[0] << " " << delta[1] << " " << delta[2] << endl;
}
void Calculate_U(atom temp_atom) {
	delta[0] = r[0] - temp_atom.r[0];
	delta[1] = r[1] - temp_atom.r[1];
	delta[2] = r[2] - temp_atom.r[2];
	d_squared = delta[1]*delta[1] + delta[2]*delta[2] + delta[0]*delta[0];
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
		string dir_name;
		float E;
		float U;
		float temp_r[3];
		float temp_v[3];
		float temp_a[3];
		int length;
	public:
		Group_Of_Atoms(int n, string filename)
		{
		dir_name = filename;
		float pos[3]{0, 0, 0};
		float zero_vec[3]{0, 0, 0};
			for (int i=0; i<n; i++)
			{
				for (int j=0; j<n; j++)
				{
					for (int k=0; k<n; k++)
					{
						pos[0] = 3*i;
						pos[1] = 3*j;
						pos[2] = 3*k;
					atom temporary_atom(pos, zero_vec, zero_vec);
					atoms.push_back(temporary_atom);
					}
				}
			}
		length = atoms.size();
		}

	
		void Group_Info()
		{
			int length = atoms.size();
			for (int i=0; i<length; i++)
			{
				atoms[i].All_Info();
			}
		}


		void One_Velocity_Verlet_Iteration()
		{
			E = 0;
			U = 0;
			vector<vector<float>> temp;
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
			for (int i=0; i<length; i++) 
			{
				t_1 = clock();
				for (int j = 0; j<length; j++)
				{
					if (i!=j)
					{
						atoms[i].Calculate_F(atoms[j]);
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
			int length = atoms.size();
			for (int i=0; i<length; i++) 
			{
				atoms[i].Update_E();
				E += atoms[i].Get_E();
				atoms[i].Set_U(0);
				for (int j = 0; j<length; j++)
				{
					if (i!=j)
					{
						atoms[i].Calculate_U(atoms[j]);
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
	Group_Of_Atoms Group_1(7, date);
	Group_1.Calculation_And_Writing_Cycle(30);
}

