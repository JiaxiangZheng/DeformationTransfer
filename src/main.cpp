#include "TriMesh.h"
#include "DeformationTransfer.h"
#include <fstream>
#include <iostream>
using namespace std;

char dataset[4][3][512] = { 
	{"./data/horse_camel/ref_horse.obj", "./data/horse_camel/ref_camel.obj", "./data/horse_camel/ref_horse_camel.cons"},
	{"face_src.obj", "face_dst.obj", "face_src_dst.cons"}, 
	{"cube_src.obj", "cube_dst.obj", "cube_src_dst.cons"},
	{"obj_src.obj", "obj_dst.obj", "obj_src_dst.cons"}
};

std::vector<bool> boundary_flag;

int main() {
	int test_index = 0;
	printf("start loading src and dst\n");
	TriMesh src(dataset[test_index][0]);
	TriMesh dst(dataset[test_index][1]);
	printf("finished loading src and dst\n");

	dst.getVertexBoundary(boundary_flag);		//���ڲ��ֶ�Ӧ�Ƚ���Ч

	vector<pair<int, int>> corres;
	ifstream ifs(dataset[test_index][2]); int cnt;
	if (!ifs) exit(-1);
	ifs >> cnt;  corres.resize(cnt);
	for (int i=0; i<cnt; ++i) ifs >> corres[i].first >> corres[i].second;

	{
		double dist = 0.0;
		for (size_t i=0; i<corres.size(); ++i) {
			Eigen::Vector3d v = src.vertex_coord[corres[i].first] - dst.vertex_coord[corres[i].second];
			dist += v.norm();
		}
		dist /= corres.size();
		std::cout << "Averge dist of the hard constraints is " << dist << std::endl;
	}
	deform_transfer(src, dst, corres);

	return 0;
}
