#ifndef ALG_CURL_H
#define ALG_CURL_H

#include <eigen/Sparse>
#include <common/meshmodel.h>
#include <common/interfaces.h>
#include <vector>
#include <math.h>
#include "alg_base.h"
using namespace std;
using namespace Eigen;
//typedef Eigen::Triplet<double> T;
//typedef Eigen::SparseMatrix<double> SparseMatrixType;

namespace gdut_curl{
	void buildRHS(MeshModel &m,std::map<int,vcg::Point3f>& map,VectorXd& b){
		// visit each face
		for(int i = 0; i < m.cm.FN(); i++) {
			CFaceO &face = m.cm.face[i];	   
			int col = face.Index();
			vcg::Point3f& kexi = map[face.Index()];
			for(int j = 0;j < 3; j++) {
				// get vertex positions
				CVertexO* v0 = face.V( (j) % 3 );
				CVertexO* v1 = face.V( (j+1) % 3 );
				CVertexO* v2 = face.V( (j+2) % 3 );

				// get vertex indices
				int row = 3*v0->Index();
				
				vcg::Point3f nabla_phi; 
				gdut_base::countPhi(v0->P(),v1->P(),v2->P(),nabla_phi);
				vcg::Point3f rhs = nabla_phi ^ kexi;
				// spin::Quaternion q(0,nabla_phi0.X()/2,nabla_phi0.Y()/2,nabla_phi0.Z()/2);
				// add contribution of this cotangent to the matrix
				b[row] += rhs.X();
				b[row+1] += rhs.Y();
				b[row+2] += rhs.Z();
			}
		}
	}
	void  buildLHS(MeshModel &m,Eigen::SparseMatrix<double> &L)
		// builds the cotan-Laplace operator
	{
		// allocate a sparse |V|x|V| matrix
		int nV = m.cm.VN();
		int nF = m.cm.FN();
		//L.reserve(VectorXi::Constant(nV,6));
		L.setZero();
		// visit each face
		for(int i = 0; i < m.cm.FN(); i++) {
			CFaceO &face = m.cm.face[i];	   
			int col = 3*face.Index();
			for(int j = 0;j < 3; j++) {
				// get vertex positions
				CVertexO* v0 = face.V( (j) % 3 );
				CVertexO* v1 = face.V( (j+1) % 3 );
				CVertexO* v2 = face.V( (j+2) % 3 );
				
				vcg::Point3f nabla_phi; 
				gdut_base::countPhi(v0->P(),v1->P(),v2->P(),nabla_phi);

				// add contribution of this cotangent to the matrix	
				int row = 3*v0->Index();
				L.coeffRef(row,col+1)-=nabla_phi.Z();
				L.coeffRef(row,col+2)+=nabla_phi.Y();

				L.coeffRef(row+1,col)+=nabla_phi.Z();
				L.coeffRef(row+1,col+2)-=nabla_phi.X();

				L.coeffRef(row+2,col)-=nabla_phi.Y();
				L.coeffRef(row+2,col+1)+=nabla_phi.X();
			}
		}
		L.makeCompressed();//压缩剩余的空间
	}



	void countDivfree(MeshModel &m,std::map<int,vcg::Point3f>& kexi,std::map<int,vcg::Point3f>& result){
		//    vector<float> x;// 解向量共|F|个分量
		int nV = m.cm.VN();
		int nF = m.cm.FN();
		VectorXd b(3*nV);
		b.setZero();
		buildRHS(m,kexi,b);

		double mean = b.mean();
		double maxb = b.maxCoeff();
		double minb = b.minCoeff();
	
		SparseMatrix<double> L(3*nV,3*nF);// 系数矩阵|V|*|F|，为（4）中左侧的部分
		L.reserve(36*nV);
		buildLHS(m,L);
		int nZ = L.nonZeros();

		Eigen::SimplicialCholesky<SparseMatrix<double>> chol(L);  // 执行A的 Cholesky分解
		//Eigen::ConjugateGradient<SparseMatrix<double>,Eigen::Upper> chol(A);  // 执行A的 Cholesky分解
		VectorXd x(3*m.cm.FN());
		x = chol.solve(b);         // 使用A的Cholesky分解来求解等号右边的向量b
		
		if(chol.info()!=Success) {
			// decomposition failed
			for(int i = 0; i < m.cm.FN();i++){
				int begin = 3*i;
				result.insert(pair<int,vcg::Point3f>(i,vcg::Point3f(x[begin],x[begin+1],x[begin+2])));
			}
			return;
		}	
	}
}

#endif