#ifndef _TRI_MESH_HH_
#define _TRI_MESH_HH_
#include <vector>
#include <cmath>
using namespace std;
#include <Eigen/Dense>

#define ERROR_LOG(ERROR_INFO)		\
	fprintf(stderr, "%s\t%s\t%d\n",	\
	ERROR_INFO, __FILE__, __LINE__)

#define DATA_09
typedef enum DisplayMode {POINT_MODE, EDGE_MODE, FACE_MODE} DisplayMode;

class TriMesh
{
public:
	struct PolyIndex
	{
		unsigned int vert_index[3];
		unsigned int norm_index[3];
	};
	vector<Eigen::Vector3d> vertex_coord;//
	vector<Eigen::Vector3d> norm_coord;
	vector<PolyIndex> polyIndex;
	vector<Eigen::Vector3d> face_norm;
	Eigen::Vector3d BoundingBox_Min, BoundingBox_Max;	//min max
	int vert_num;
	int poly_num;

	TriMesh() : vert_num(0), poly_num(0) {}
	TriMesh(const char* filename);
	void updateNorm();
	void render(DisplayMode mode) ;
	void prepareFaceNeighbours(std::vector<std::vector<int>>& neighbours) const;
	void getVertexBoundary(std::vector<bool>& v_flags) const;
	void saveOBJ(const char* filename) const;
	void getBoundingBox(Eigen::Vector3d& Min, Eigen::Vector3d& Max);
private:
	void readPLY(const char* filename);
	void readOBJ(const char* filename);
};

#endif/*_TRI_MESH_HH_*/



