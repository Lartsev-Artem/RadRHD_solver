#include "reader_txt.h"

int ReadSphereDirectionDecart(const std::string name_file_sphere_direction, std::vector<Vector3>& directions_all) {

	std::ifstream ifile;

	OPEN_FSTREAM(ifile, name_file_sphere_direction.c_str());

	int N = 0;
	ifile >> N;
	directions_all.resize(N);

	Type buf_s;
	int i = 0;

	for (int i = 0; i < N; i++) {
		ifile >> buf_s;
		ifile >> directions_all[i][0] >> directions_all[i][1] >> directions_all[i][2];
	}
	ifile >> buf_s;
	ifile.close();
	return 0;
}

size_t ReadSphereDirectionDecartToSpherical(const std::string& name_file_sphere_direction, grid_directions_t& grid_direction)
{

	std::ifstream ifile;
	OPEN_FSTREAM(ifile, name_file_sphere_direction.c_str());

	int N = 0;
	ifile >> N;

	grid_direction.size = N;
	grid_direction.directions.resize(N);

	for (int i = 0; i < N; i++) {
		ifile >> grid_direction.directions[i].area;
		ifile >> grid_direction.directions[i].dir[0] 
			  >> grid_direction.directions[i].dir[1] 
			  >> grid_direction.directions[i].dir[2];
	}
	ifile >> grid_direction.full_area;
	ifile.close();

	return 0;
}
