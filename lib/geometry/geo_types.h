#ifndef GEO_TYPES
#define GEO_TYPES

#if defined LINUX
#include <eigen3/Eigen/Dense>
#elif defined WINDOWS
#include "C:/DEV/vcpkg-master/installed/x64-windows/include/eigen3/Eigen/Dense"
#elif defined CLASTER
#include </nethome/student/FS18/FS2-x1/Lartsev/Eigen/Dense>
#endif

#include "global_def.h"

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::VectorXd VectorX;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::Vector4d Vector4;
typedef Eigen::Matrix4d Matrix4;
typedef Eigen::MatrixXd MatrixX;
typedef double Type;
typedef uint8_t ShortId;


///\todo: d в новый файл  с классами
struct Normals {
	std::vector<Vector3> n;
	Normals(const int size = CELL_SIZE) { n.resize(size); }
	~Normals() { n.clear(); }
};

struct Face {
	Vector3 A;
	Vector3 B;
	Vector3 C;
	Face& operator=(const Face& face) 
	{
		A = face.A;
		B = face.B;
		C = face.C;
		return *this;
	}
};

struct FaceCell {
	int face_id;
	Face face;
	FaceCell(const int id = 0, const Face& face_init = Face()) 
	: face_id(id), face(face_init) {}
};

struct direction_s
{
	Vector3 dir;
	Type area;
};

struct grid_directions_t
{
	int size;
	std::vector<direction_s> directions;
	Type full_area;
};

struct BasePointTetra //узлы интерполяции всех тетраэдров // В перспективе можно уйти к граням
{
	Vector3 x[CELL_SIZE][NUMBER_OF_MEASUREMENTS];

	Vector3 operator()(const int i, const int j) const
	{
		return x[i][j];
	}
};
struct cell_local // для каждой ячейки и каждого направления
{
	Type s;    //расстояние x0-x
	Vector2 x0; //локальная координата входного узла для интерполяции
	ShortId in_face_id; //id выходной грани

	friend std::ostream& operator<< (std::ostream& out, const cell_local& point);
};

struct grid_t
{
	struct elem_t
	{
		int a;
	};

	struct face_t
	{
		int a;
	};
	

	int size;	
	std::vector<elem_t> cells;
	std::vector<face_t> faces;

	/// \todo all config!
};

#endif //GEO_TYPES