#ifndef ALG_CURLFREE_H
#define ALG_CURLFREE_H

#include <eigen/Sparse>
#include <common/meshmodel.h>
#include <common/interfaces.h>
#include <vector>
#include <math.h>
using namespace std;
using namespace Eigen;
//typedef Eigen::Triplet<double> T;
//typedef Eigen::SparseMatrix<double> SparseMatrixType;

namespace gdut_curl_free{
	
	void buildB(MeshModel &m,std::map<int,vcg::Point3f>& map,VectorXd& b){
		// visit each face
		for(int i = 0; i < m.cm.FN(); i++) {
			CFaceO &face = m.cm.face[i];	   
			vcg::Point3f& kexi = map[face.Index()];
			float area = gdut_base::countArea(face);
			vcg::Point3f normalized = ::NormalizedNormal(face);
			for(int j = 0;j < 3; j++) {
				// get vertex positions
				CVertexO* v0 = face.V( (j) % 3 );
				CVertexO* v1 = face.V( (j+1) % 3 );
				CVertexO* v2 = face.V( (j+2) % 3 );
				int row = v0->Index();		

				vcg::Point3f nabla_phi0;
				gdut_base::countPhi(v0->P(),v1->P(),v2->P(),nabla_phi0);

				b[row] += (nabla_phi0*kexi)*area;
			}
		}
	}
	//void  buildLaplacian(MeshModel &m,Eigen::SparseMatrix<double> &L)
	//{
	//	int nV = m.cm.VN();
	//	int nF = m.cm.FN();
	//	L.reserve(VectorXi::Constant(nV,6));
	//	L.setZero();
	//	// visit each face
	//	for(int i = 0; i < nF; i++) {
	//		CFaceO &face = m.cm.face[i];	   
	//		int col = face.Index();
	//		for(int j = 0;j < 3; j++) {
	//			// get vertex positions
	//			CVertexO* v0 = face.V( (j) % 3 );
	//			CVertexO* v1 = face.V( (j+1) % 3 );
	//			CVertexO* v2 = face.V( (j+2) % 3 );

	//			// get vertex indices
	//			int k1 = v1->Index();
	//			int k2 = v2->Index();

	//			// compute cotangent of the angle at the current vertex
	//			// (equal to cosine over sine, which equals the dot
	//			// product over the norm of the cross product)       
	//			vcg::Point3f u1 = v1->P() - v0->P();
	//			vcg::Point3f u2 = v2->P() - v0->P();
	//			double cotAlpha = (u1*u2)/(u1^u2).Norm();

	//			// add contribution of this cotangent to the matrix
	//			/* L( k1, k2 ) -= cotAlpha / 2.;
	//			L( k2, k1 ) -= cotAlpha / 2.;
	//			L( k1, k1 ) += cotAlpha / 2.;
	//			L( k2, k2 ) += cotAlpha / 2.;*/

	//			L.coeffRef(k1,k2)+=(-cotAlpha / 2.);
	//			L.coeffRef(k2,k1)+=(-cotAlpha / 2.);
	//			L.coeffRef(k1,k1)+=(cotAlpha / 2.);
	//			L.coeffRef(k2,k2)+=(-cotAlpha / 2.);
	//		}
	//	}
	//	L.makeCompressed();//压缩剩余的空间
	//}

	void  buildLaplacian(MeshModel &m,Eigen::SparseMatrix<double> &L)
	{
		int nV = m.cm.VN();
		int nF = m.cm.FN();
		L.reserve(VectorXi::Constant(nV,6));
		L.setZero();
		// visit each face
		for(int i = 0; i < nF; i++) {
			CFaceO &face = m.cm.face[i];	   
			float area = gdut_base::countArea(face);

				// get vertex positions
				CVertexO* v0 = face.V(0);
				CVertexO* v1 = face.V(1);
				CVertexO* v2 = face.V(2);

				// get vertex indices
				int k0 = v0->Index();
				int k1 = v1->Index();
				int k2 = v2->Index();

				vcg::Point3f nabla_phi0; 
				gdut_base::countPhi(v0->P(),v1->P(),v2->P(),nabla_phi0);
				vcg::Point3f nabla_phi1; 
				gdut_base::countPhi(v1->P(),v2->P(),v0->P(),nabla_phi1);
				vcg::Point3f nabla_phi2; 
				gdut_base::countPhi(v2->P(),v0->P(),v1->P(),nabla_phi2);

				float n01 = nabla_phi0*nabla_phi1;
				float n12 = nabla_phi1*nabla_phi2;
				float n02 = nabla_phi0*nabla_phi2;

				L.coeffRef(k0,k0)+=nabla_phi0*nabla_phi0*area;
				L.coeffRef(k0,k1)+=n01*area;
				L.coeffRef(k0,k2)+=n02*area;
				L.coeffRef(k1,k0)+=n01*area;
				L.coeffRef(k1,k1)+=nabla_phi1*nabla_phi1*area;
				L.coeffRef(k1,k2)+=n12*area;
				L.coeffRef(k2,k0)+=n02*area;
				L.coeffRef(k2,k1)+=n12*area;
				L.coeffRef(k2,k2)+=nabla_phi2*nabla_phi2*area;
			
		}
		L.makeCompressed();//压缩剩余的空间
	}

	void count_nabla_U(CFaceO &f,VectorXd& x,vcg::Point3f& u){
		CVertexO* v0 = f.V(0);
				CVertexO* v1 = f.V(1);
				CVertexO* v2 = f.V(2);
				double s0=x[v0->Index()];
				double s1=x[v1->Index()];
				double s2=x[v2->Index()];
				vcg::Point3f nablaD;
				gdut_base::countNablaOfFace(f,s0,s1,s2,u);
	}

	void countCurlfree(MeshModel &m,std::map<int,vcg::Point3f>& kexi,VectorXd& x){

		int nV = m.cm.VN();
		VectorXd b(nV);
		b.setZero();
		buildB(m,kexi,b);

		double mean = b.mean();
		double maxb = b.maxCoeff();
		double minb = b.minCoeff();
		
		SparseMatrix<double> L(nV,nV);// 系数矩阵|V|*|V|
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