#include <vector>
#include <set> 
#include <map>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "DeformationTransfer.h"
#include "TriMesh.h"
using namespace Eigen;
static void _solve_linear_system(TriMesh& src, const TriMesh& dst, 
    const std::vector<std::vector<int>>& tri_neighbours, 
    const std::vector<std::pair<int, int> >& corres,
    double weights[3]
);
static void _build_kdtree(const TriMesh& dataset_mesh);
static void _release_kdtree();
static void _setup_kd_correspondence(
    const TriMesh& src_mesh, const TriMesh& dst_mesh,
    std::vector<std::pair<int,int>>& soft_correspondence
);

struct edge_t {
    unsigned int _v[2];
    edge_t(unsigned int v1, unsigned int v2) {
        _v[0] = v1 > v2 ? v2 : v1;
        _v[1] = v1 > v2 ? v1 : v2;
    }
};
inline bool 
operator< (const edge_t& ref1, const edge_t& ref2) {
        return (ref1._v[0] < ref2._v[0] ? true : (ref1._v[0] == ref2._v[0] && ref1._v[1] < ref2._v[1]));
}

static void 
_prepare_neighbours(const TriMesh& mesh, std::vector<std::vector<int>>& tri_neighbours) {
	int tri_num = mesh.poly_num; 
    std::multimap<edge_t, int> edges;
    std::set<edge_t> edge_sets;
    for (int i=0; i<tri_num; ++i) {
        const unsigned int* v_index = mesh.polyIndex[i].vert_index;
        for (int j=0; j<3; ++j) {
            edge_t edge(v_index[j], v_index[(j+1)%3]);
            edges.insert(std::make_pair(edge, i));
            edge_sets.insert(edge);
        }
    }
    
    tri_neighbours.resize(tri_num);
    for (auto it = edge_sets.begin(); it != edge_sets.end(); ++it) {
        auto lb = edges.lower_bound(*it);
        auto ub = edges.upper_bound(*it);

        std::vector<int> face_indexs;
        for (; lb != edges.end() && lb != ub; ++lb)
            face_indexs.push_back(lb->second);
        for (auto it=face_indexs.begin(); it != face_indexs.end(); ++it) {
            for (auto it_ = face_indexs.begin(); it_ != face_indexs.end(); ++it_) {
                if (it == it_) continue;
                tri_neighbours[*it].push_back(*it_);
            }
        }
    }

    return;
}

void 
deform_transfer(TriMesh& src, const TriMesh& dst, const std::vector<std::pair<int, int> >& corres
) {
    std::vector<std::vector<int>> tri_neighbours;
    _prepare_neighbours(src, tri_neighbours);
    _build_kdtree(dst);

    double weights[4] = {1.0, 10000, 0.1, 0.0};  //smooth, regular, identity weights, soft_constraints(nearest point)
    
	int cnt = 0;    //change the step vale and upper bound of w[3] according to the data.
	for (; weights[3] <= 200.0; weights[3]+=20.0, ++cnt) {
		_solve_linear_system(src, dst, tri_neighbours, corres, weights);
		char filename[512]; sprintf(filename, "src_to_dst_%d.obj", cnt);
		src.saveOBJ(filename);
	}

	_release_kdtree();
    return;
}

static void 
_solve_linear_system(TriMesh& src, const TriMesh& dst, 
    const std::vector<std::vector<int>>& tri_neighbours, 
    const std::vector<std::pair<int, int> >& corres,
    double weights[4]
) {
    int rows = 0; int cols = src.vert_num + src.poly_num;
    assert (tri_neighbours.size() == src.poly_num);
    std::vector<std::pair<int, int>> soft_corres;
    _setup_kd_correspondence(src, dst, soft_corres);

    for (int i=0; i<src.poly_num; ++i)
        rows += 3*tri_neighbours[i].size();
    rows += src.poly_num*3;
    rows += corres.size();
    rows += soft_corres.size();
    SparseMatrix<double> A(rows, cols);  
	A.reserve(Eigen::VectorXi::Constant(cols, 200));
    std::vector<VectorXd> Y(3, VectorXd(rows));
    Y[0].setZero();Y[1].setZero();Y[2].setZero();

	src.updateNorm();
    //precompute Q_hat^-1 in Q_hat_inverse   [n v2-v1 v3-v1]^-1
    assert (!src.face_norm.empty());
    std::vector<Matrix3d> Q_hat_inverse(src.poly_num);
	Eigen::Matrix3d _inverse;
    for (int i=0; i<src.poly_num; ++i) {
        unsigned int* index_i = src.polyIndex[i].vert_index;
        Vector3d v[3];
		v[0] = src.face_norm[i];
        v[1] = src.vertex_coord[index_i[1]] - src.vertex_coord[index_i[0]];//v2-v1
        v[2] = src.vertex_coord[index_i[2]] - src.vertex_coord[index_i[0]];//v3-v1
        for (int k=0; k<3; ++k)
            for (int j=0; j<3; ++j)
                Q_hat_inverse[i](k, j) = v[j][k];
		//std::cout << Q_hat_inverse[i] << std::endl;
		_inverse = Q_hat_inverse[i].inverse();
		Q_hat_inverse[i] = _inverse;
    }
	int energy_size[4] = {0, 0, 0, 0};

    //start establishing the large linear sparse system
    double weight_smooth = weights[0]; int row = 0;
    for (int i=0; i<src.poly_num; ++i) {
        Eigen::Matrix3d& Q_i_hat = Q_hat_inverse[i];
        unsigned int* index_i = src.polyIndex[i].vert_index;
        for (size_t _j=0; _j<tri_neighbours[i].size(); ++_j) {
            int j = tri_neighbours[i][_j];
            Eigen::Matrix3d& Q_j_hat = Q_hat_inverse[j];
            unsigned int* index_j = src.polyIndex[j].vert_index;
			for (int k=0; k<3; ++k) {
				for (int p=0; p<3; ++p)
					A.coeffRef(row+k, src.poly_num+index_i[p]) = 0.0;
			}
            for (int k=0; k<3; ++k) {//for each row
                A.coeffRef(row + k, i) = weight_smooth*Q_i_hat(0, k);   //n
                A.coeffRef(row+k, src.poly_num + index_i[0]) += 
                    -weight_smooth*(Q_i_hat(1, k) + Q_i_hat(2, k));
                A.coeffRef(row+k, src.poly_num + index_i[1]) += 
                    weight_smooth*Q_i_hat(1, k);
                A.coeffRef(row+k, src.poly_num + index_i[2]) += 
                    weight_smooth*Q_i_hat(2, k);

                A.coeffRef(row+k, j) = -weight_smooth*Q_j_hat(0, k);   //n
                A.coeffRef(row+k, src.poly_num + index_j[0]) += 
                    weight_smooth*(Q_j_hat(1, k) + Q_j_hat(2, k));
                A.coeffRef(row+k, src.poly_num + index_j[1]) += 
                    -weight_smooth*Q_j_hat(1, k);
                A.coeffRef(row+k, src.poly_num + index_j[2]) += 
                    -weight_smooth*Q_j_hat(2, k);
            }
			row += 3;
        }
    }
	energy_size[0] = row;

    double weight_regular = weights[1];
    for (int i=0; i<corres.size(); ++i, ++row) {
		A.coeffRef(row, src.poly_num + corres[i].first) = weight_regular;
        for (int j=0; j<3; ++j)
            Y[j](row) = weight_regular*dst.vertex_coord[corres[i].second][j];
    }
	energy_size[1] = row;

	double weight_identity = weights[2];
    for (int i=0; i<src.poly_num; ++i) {
        Eigen::Matrix3d Q_i_hat = Q_hat_inverse[i];
        unsigned int* index_i = src.polyIndex[i].vert_index;
        Y[0](row) = weight_identity;    Y[0](row+1) = 0.0;              Y[0](row+2) = 0.0; 
        Y[1](row) = 0.0;                Y[1](row+1) = weight_identity;  Y[1](row+2) = 0.0; 
        Y[2](row) = 0.0;                Y[2](row+1) = 0.0;              Y[2](row+2) = weight_identity; 
        for (int k=0; k<3; ++k, ++row) {
            A.coeffRef(row, i) = weight_identity*Q_i_hat(0, k);   //n
            A.coeffRef(row, src.poly_num + index_i[0]) = 
                -weight_identity*(Q_i_hat(1, k) + Q_i_hat(2, k));
            A.coeffRef(row, src.poly_num + index_i[1]) = 
                weight_identity*Q_i_hat(1, k);
            A.coeffRef(row, src.poly_num + index_i[2]) = 
                weight_identity*Q_i_hat(2, k);
        }
    }
	energy_size[2] = row;

    double weight_soft_constraint = weights[3];
    for (int i=0; i<soft_corres.size(); ++i, ++row) {
		A.coeffRef(row, src.poly_num + soft_corres[i].first) = weight_soft_constraint;
        for (int j=0; j<3; ++j)
            Y[j](row) = weight_soft_constraint*dst.vertex_coord[soft_corres[i].second][j];
    }
	energy_size[3] = row;

   	//start solving the least-square problem
	fprintf(stdout, "finished filling matrix\n");
	Eigen::SparseMatrix<double> At = A.transpose();
	Eigen::SparseMatrix<double> AtA = At*A;

	Eigen::SimplicialCholesky<SparseMatrix<double>> solver;
	solver.compute(AtA);
	if (solver.info() != Eigen::Success) {
		fprintf(stdout, "unable to defactorize AtA\n");
		exit(-1);
	}

    VectorXd X[3];
    for (int i=0; i<3; ++i) {
        VectorXd AtY = At*Y[i];
        X[i] = solver.solve(AtY);
		Eigen::VectorXd Energy = A*X[i] - Y[i];
		Eigen::VectorXd smoothEnergy = Energy.head(energy_size[0]);
		Eigen::VectorXd hardRegularEnergy = Energy.segment(energy_size[0], energy_size[1]-energy_size[0]);
		Eigen::VectorXd identityEnergy = Energy.segment(energy_size[1], energy_size[2]-energy_size[1]);
		Eigen::VectorXd softRegularEnergy = Energy.tail(energy_size[3]-energy_size[2]);
		fprintf(stdout, "\t%lf = %lf + %lf + %lf + %lf\n", 
			Energy.dot(Energy), smoothEnergy.dot(smoothEnergy), hardRegularEnergy.dot(hardRegularEnergy), 
			identityEnergy.dot(identityEnergy), softRegularEnergy.dot(softRegularEnergy));
    }
    
    //fill data back to src
    for (int i=0; i<src.poly_num; ++i)
        for (int d=0; d<3; ++d)
            src.face_norm[i][d] = X[d](i);
    for (int i=0; i<src.vert_num; ++i) 
        for (int d=0; d<3; ++d)
            src.vertex_coord[i][d] = X[d](src.poly_num+i);
    
    return;
}

#include <flann/flann.hpp>
static double* dataset;
static flann::Index<flann::L2<double>> *kd_flann_index;

//the threshold should be adjust according to the input data or you can just use
//the hard constraints's average error to compute the threshold
#define THRESHOLD_DIST 0.1     //less
#define THRESHOLD_NORM 0.7     //greater
static void 
_build_kdtree(const TriMesh& dataset_mesh) {
    int vert_num = dataset_mesh.vert_num;
    dataset = new double[3*vert_num];
    {
        int _index = 0;
        for (auto it = dataset_mesh.vertex_coord.begin();
            it != dataset_mesh.vertex_coord.end(); ++it) {
                dataset[_index] = (*it)[0];
                dataset[_index+1] = (*it)[1];
                dataset[_index+2] = (*it)[2];
                _index += 3;
            }
    }
	flann::Matrix<double> flann_dataset(dataset, vert_num, 3);
    kd_flann_index = new flann::Index<flann::L2<double> >(flann_dataset, flann::KDTreeIndexParams(1)); 
    kd_flann_index->buildIndex();
    return;
}

static void
_setup_kd_correspondence(const TriMesh& src_mesh, const TriMesh& dataset_mesh, 
std::vector<std::pair<int,int>>& soft_correspondence) {
    const int dim = 3; const int knn = 1;
    //pre allocate 
    flann::Matrix<double> query(new double[dim], 1, dim);
   	flann::Matrix<int> indices(new int[query.rows*knn], query.rows, knn);
	flann::Matrix<double> dists(new double[query.rows*knn], query.rows, knn); 

    soft_correspondence.clear();
    for (int i=0; i<src_mesh.vert_num; ++i) {
        for (int k=0; k<3; ++k) query[0][k] = src_mesh.vertex_coord[i][k];
        kd_flann_index->knnSearch(query, indices, dists, 1, 
			flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
        if (dists[0][0] < THRESHOLD_DIST && 
			src_mesh.norm_coord[i].dot(dataset_mesh.norm_coord[indices[0][0]]) >= THRESHOLD_NORM) {
            soft_correspondence.push_back(std::make_pair(i, indices[0][0]));
        }
    }

    delete []query.ptr();
    delete []indices.ptr();
    delete []dists.ptr();
    return;
}

static void
_release_kdtree() {
    delete []dataset;
	delete kd_flann_index;
	kd_flann_index = NULL;
}
