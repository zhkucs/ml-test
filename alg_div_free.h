#ifndef ALG_DIV_FREE_H
#define ALG_DIV_FREE_H

#include <eigen/Dense>
#include <common/meshmodel.h>
#include <common/interfaces.h>
#include <unordered_map>
#include <vector>
#include <math.h>
#include "alg_base.h"
using namespace std;
using namespace Eigen;
//typedef Eigen::Triplet<double> T;
//typedef Eigen::SparseMatrix<double> SparseMatrixType;

namespace gdut_div_free{
	void buildRHS(MeshModel &m,std::map<int,vcg::Point3f>& map,VectorXd& b){
		// visit each face
		for(int i = 0; i < m.cm.FN(); i++) {
			CFaceO &face = m.cm.face[i];	   
			int col = face.Index();
			float area = gdut_base::countArea(face);
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
				vcg::Point3f rhs = (nabla_phi ^ kexi)*area;
				// spin::Quaternion q(0,nabla_phi0.X()/2,nabla_phi0.Y()/2,nabla_phi0.Z()/2);
				// add contribution of this cotangent to the matrix
				b[row] += rhs.X();
				b[row+1] += rhs.Y();
				b[row+2] += rhs.Z();
			}
		}
	}

	void buildN2N(vcg::Point3f& a,vcg::Point3f& b,Matrix3d& A){
		float a1 = a.X();
		float a2 = a.Y();
		float a3 = a.Z();
		float b1 = b.X();
		float b2 = b.Y();
		float b3 = b.Z();

		A.setZero();
		A << -(a2*b2+a3*b3), a2*b1 ,       a3*b1,
			 -a1*b2,         a1*b1+a3*b3,  -a3*b2,
			a1*b3,           a2*b3,        -(a1*b1+a2*b2);
	}

	void set(Eigen::SparseMatrix<double> &L,int vi,int vj,Matrix3d& m){
		int row=vi*3;
		int col = vj*3;

		for(int i = 0; i < 3;i++){
			for(int j = 0; j < 3;j++){
				L.coeffRef(row+i,col+j) += m.coeff(i,j);
			}
		}
	}

	void  buildLHS(MeshModel &m,Eigen::SparseMatrix<double> &L)// builds the cotan-Laplace operator
	{
		// allocate a sparse |V|x|V| matrix
		int nV = m.cm.VN();
		int nF = m.cm.FN();
		// 计算各v点在面f内的梯度
		std::unordered_map<std::pair<int, int>, vcg::Point3f, gdut_base::pairhash>  nabla_v_face;		
		for(int i = 0; i < nF; i++) {
			CFaceO &face = m.cm.face[i];	
			int f_idx = face.Index();			
			// get vertex positions
			CVertexO* v0 = face.V(0);
			CVertexO* v1 = face.V(1);
			CVertexO* v2 = face.V(2);

			int i0 = v0->Index();
			int i1 = v1->Index();
			int i2 = v2->Index();

			vcg::Point3f nabla_phi0; 
			gdut_base::countPhi(v0->P(),v1->P(),v2->P(),nabla_phi0);

			vcg::Point3f nabla_phi1; 
			gdut_base::countPhi(v1->P(),v2->P(),v0->P(),nabla_phi1);

			vcg::Point3f nabla_phi2; 
			gdut_base::countPhi(v2->P(),v0->P(),v1->P(),nabla_phi2);

			float area = gdut_base::countArea(face);
			nabla_v_face[std::make_pair(i0,f_idx)]=nabla_phi0*area;
			nabla_v_face[std::make_pair(i1,f_idx)]=nabla_phi1*area;
			nabla_v_face[std::make_pair(i2,f_idx)]=nabla_phi2*area;
		}

		L.setZero();
		for(int i = 0; i < nF; i++) {
			CFaceO &face = m.cm.face[i];	   
			int f_idx = face.Index();
			CVertexO* v0 = face.V(0);
			CVertexO* v1 = face.V(1);
			CVertexO* v2 = face.V(2);
			int v_idx[3];
			v_idx[0] = v0->Index();
			v_idx[1] = v1->Index();
			v_idx[2] = v2->Index();

			vcg::Point3f na[3];
			na[0] = nabla_v_face[std::make_pair(v_idx[0],f_idx)];
			na[1] = nabla_v_face[std::make_pair(v_idx[1],f_idx)];
			na[2] = nabla_v_face[std::make_pair(v_idx[2],f_idx)];

			Matrix3d m[3][3];
			for(int i = 0; i < 3;i++){
				vcg::Point3f ni = na[i];
				for(int j = 0; j < 3;j++){
					vcg::Point3f nj = na[j];	
					buildN2N(ni,nj,m[i][j]);
					set(L,v_idx[i],v_idx[j],m[i][j]);

				}
			}
		}	
	//  L.makeCompressed();//压缩剩余的空间
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

	SparseMatrix<double> L(3*nV,3*nV);// 系数矩阵|3V|*|3V|，为（4）中左侧的部分
	L.reserve(36*nV);
	int nZ1 = L.nonZeros();
	buildLHS(m,L);
	int nZ = L.nonZeros();

	Eigen::SimplicialCholesky<SparseMatrix<double>> chol(L);  // 执行A的 Cholesky分解
	//Eigen::ConjugateGradient<SparseMatrix<double>,Eigen::Upper> chol(A);  // 执行A的 Cholesky分解
	VectorXd x(3*m.cm.VN());
	x = chol.solve(b);         // 使用A的Cholesky分解来求解等号右边的向量b

	double meanX = x.mean();
	double maxbX = x.maxCoeff();
	double minbX = x.minCoeff();

	if(chol.info()==Success) {
		for(int i = 0; i < nF;i++){
			CFaceO face = m.cm.face[i];
			CVertexO* v0 = face.V(0);
			CVertexO* v1 = face.V(1);
			CVertexO* v2 = face.V(2);

			int i0 = v0->Index()*3;
			int i1 = v1->Index()*3;
			int i2 = v2->Index()*3;

			int begin = 3*i;
			vcg::Point3f p0(x[i0],x[i0+1],x[i0+2]);
			vcg::Point3f p1(x[i1],x[i1+1],x[i1+2]);
			vcg::Point3f p2(x[i2],x[i2+1],x[i2+2]);

			vcg::Point3f nabla_phi0; 
			gdut_base::countPhi(v0->P(),v1->P(),v2->P(),nabla_phi0);

			vcg::Point3f nabla_phi1; 
			gdut_base::countPhi(v1->P(),v2->P(),v0->P(),nabla_phi1);

			vcg::Point3f nabla_phi2; 
			gdut_base::countPhi(v2->P(),v0->P(),v1->P(),nabla_phi2);

			vcg::Point3f v_f = nabla_phi0 ^ p0 + nabla_phi1 ^ p1 + nabla_phi2 ^ p2;// 各顶点与基函数梯度叉积后求和，为三角形面片内
			result.insert(pair<int,vcg::Point3f>(i,v_f));
		}
		return;
	}	
}
}

#endif