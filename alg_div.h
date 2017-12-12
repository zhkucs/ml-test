#ifndef ALG_EIGEN_H
#define ALG_EIGEN_H

#include <eigen/Sparse>
#include <common/meshmodel.h>
#include <common/interfaces.h>
#include <vector>
#include <math.h>
using namespace std;
using namespace Eigen;
//typedef Eigen::Triplet<double> T;
//typedef Eigen::SparseMatrix<double> SparseMatrixType;

namespace gdut{
	
	void buildB(MeshModel &m,std::map<int,vcg::Point3f>& map,VectorXd& b){
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
				int row = v0->Index();		

				// compute cotangent of the angle at the current vertex
				// (equal to cosine over sine, which equals the dot
				// product over the norm of the cross product)
				vcg::Point3f u12 = v2->P() - v1->P();
				vcg::Point3f u01 = v1->P() - v0->P();
				vcg::Point3f normal_face=(u01^u12).Normalize();

				vcg::Point3f nabla_phi0 = u12^normal_face;
				// spin::Quaternion q(0,nabla_phi0.X()/2,nabla_phi0.Y()/2,nabla_phi0.Z()/2);
				// add contribution of this cotangent to the matrix
				b[row] += nabla_phi0*kexi;
			}
		}
	}
	void  buildLaplacian(MeshModel &m,Eigen::SparseMatrix<double> &L)
		// builds the cotan-Laplace operator
	{
		// allocate a sparse |V|x|V| matrix
		int nV = m.cm.VN();
		int nF = m.cm.FN();
		L.reserve(VectorXi::Constant(nV,6));
		L.setZero();
		// visit each face
		for(int i = 0; i < m.cm.FN(); i++) {
			CFaceO &face = m.cm.face[i];	   
			int col = face.Index();
			for(int j = 0;j < 3; j++) {
				// get vertex positions
				CVertexO* v0 = face.V( (j) % 3 );
				CVertexO* v1 = face.V( (j+1) % 3 );
				CVertexO* v2 = face.V( (j+2) % 3 );

				// get vertex indices
				int k1 = v1->Index();
				int k2 = v2->Index();

				// compute cotangent of the angle at the current vertex
				// (equal to cosine over sine, which equals the dot
				// product over the norm of the cross product)       
				vcg::Point3f u1 = v1->P() - v0->P();
				vcg::Point3f u2 = v2->P() - v0->P();
				double cotAlpha = (u1*u2)/(u1^u2).Norm();

				// add contribution of this cotangent to the matrix
				/* L( k1, k2 ) -= cotAlpha / 2.;
				L( k2, k1 ) -= cotAlpha / 2.;
				L( k1, k1 ) += cotAlpha / 2.;
				L( k2, k2 ) += cotAlpha / 2.;*/

				L.coeffRef(k1,k2)+=(-cotAlpha / 2.);
				L.coeffRef(k2,k1)+=(-cotAlpha / 2.);
				L.coeffRef(k1,k1)+=(cotAlpha / 2.);
				L.coeffRef(k2,k2)+=(-cotAlpha / 2.);
			}
		}
		L.makeCompressed();//压缩剩余的空间
	}



	void countCurlfree(MeshModel &m,std::map<int,vcg::Point3f>& kexi,VectorXd& x){
		//    vector<float> x;// 解向量共|F|个分量
		VectorXd b(m.cm.VN());
		b.setZero();
		buildB(m,kexi,b);

		double mean = b.mean();
		double maxb = b.maxCoeff();
		double minb = b.minCoeff();

		int nV = m.cm.VN();
		SparseMatrix<double> L(m.cm.VN(),m.cm.VN());// 系数矩阵|V|*|F|，为（4）中左侧的部分
		L.reserve(6*nV);
		buildLaplacian(m,L);
		int nZ = L.nonZeros();

		Eigen::SimplicialCholesky<SparseMatrix<double>> chol(L);  // 执行A的 Cholesky分解
		//Eigen::ConjugateGradient<SparseMatrix<double>,Eigen::Upper> chol(A);  // 执行A的 Cholesky分解
		x = chol.solve(b);         // 使用A的Cholesky分解来求解等号右边的向量b
		
		if(chol.info()!=Success) {
			// decomposition failed
			std::cout  << x <<std::endl;
			double mm = x.mean();
			double max = x.maxCoeff();
			double min = x.minCoeff();
			return;
		}	
	}
}

#endif