#include "TriMesh.h"
#include "DeformationTransfer.h"
#include <fstream>
#include <iostream>
using namespace std;

char dataset[6][3][512] = { 
	{"ref_horse.obj", "ref_camel.obj", "ref_horse_camel.cons"},
	{"face_src.obj", "face_dst.obj", "face_src_dst.cons"}, 
	{"cube_src.obj", "cube_dst.obj", "cube_src_dst.cons"},
	{"obj_src.obj", "obj_dst.obj", "obj_src_dst.cons"},
	{	"E:/data/ScapeOriginal/build/data_ming/pview_front_upsampling.obj", 
		"E:/data/ScapeOriginal/build/data_ming/Register_front.obj",
		"E:/data/ScapeOriginal/build/data_ming/frontPair_upsampling.txt"
	}, 
	{
		"F:/temp/pview_front_upsampling.obj",
		"F:/temp/merged_result.obj",
		"F:/temp/corres.cons"
	}
};

namespace {
	struct edge_t {
		unsigned int v1, v2;
		edge_t(unsigned int _v1, unsigned int _v2) {
			assert (_v1 != _v2);
			v1 = _v1 < _v2 ? _v1 : _v2;
			v2 = _v1 < _v2 ? _v2 : _v1;
		}
	};
	bool operator < (const edge_t& ref1, const edge_t& ref2) {
		return ref1.v1 < ref2.v1 ? true : (ref1.v1 == ref2.v1 && ref1.v2 < ref2.v2);
	}

}//end namespace
#include <map>
static void mesh_boundary_detect(const TriMesh& mesh, std::vector<bool>& boundary_v_flags) {
	int nv = mesh.vert_num, nf = mesh.poly_num;
	std::map<edge_t, int> edges;
	for (int i=0; i<nf; ++i) {
		edges[edge_t(mesh.polyIndex[i].vert_index[0], mesh.polyIndex[i].vert_index[1])]++;
		edges[edge_t(mesh.polyIndex[i].vert_index[1], mesh.polyIndex[i].vert_index[2])]++;
		edges[edge_t(mesh.polyIndex[i].vert_index[2], mesh.polyIndex[i].vert_index[0])]++;
	}
	boundary_v_flags.clear(); boundary_v_flags.resize(nv, false);
	for (auto it=edges.begin(); it != edges.end(); ++it) {
		if (it->second != 2) {
			boundary_v_flags[(it->first).v1] = true;
			boundary_v_flags[(it->first).v2] = true;
		}
	}

	return;
}
std::vector<bool> boundary_flag;

int main() {
	int test_index = 5;
 	mesh_boundary_detect(dataset[test_index][1], boundary_flag);

	printf("start loading src and dst\n");
	TriMesh src(dataset[test_index][0]);
	TriMesh dst(dataset[test_index][1]);
	printf("finished loading src and dst\n");

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
