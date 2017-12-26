#ifndef ALG_CURL_H
#define ALG_CURL_H

#include <eigen/Sparse>
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

namespace gdut_curl{
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
	void  buildLHS(MeshModel &m,Eigen::SparseMatrix<double> &L)// builds the cotan-Laplace operator
	{
		// allocate a sparse |V|x|V| matrix
		int nV = m.cm.VN();
		std::unordered_map<std::pair<int, int>, vcg::Point3f, gdut_base::pairhash>  LL;
		//L.reserve(VectorXi::Constant(nV,6));
		// visit each face set zero
		//for(int i = 0; i < m.cm.FN(); i++) {
		//	CFaceO &face = m.cm.face[i];		
		//	// get vertex positions
		//	CVertexO* v0 = face.V(0);
		//	CVertexO* v1 = face.V(1);
		//	CVertexO* v2 = face.V(2);
		//	int i0 = v0->Index();
		//	int i1 = v1->Index();
		//	int i2 = v2->Index();

		//	LL[std::make_pair(i0,i1)].SetZero();
		//	LL[std::make_pair(i1,i0)].SetZero();
		//	LL[std::make_pair(i1,i2)].SetZero();
		//	LL[std::make_pair(i2,i1)].SetZero();
		//	LL[std::make_pair(i2,i0)].SetZero();
		//	LL[std::make_pair(i0,i2)].SetZero();

		//	LL[std::make_pair(i0,i0)].SetZero();
		//	LL[std::make_pair(i1,i1)].SetZero();
		//	LL[std::make_pair(i2,i2)].SetZero();
		//}
		for(int i = 0; i < m.cm.FN(); i++) {
			CFaceO &face = m.cm.face[i];	   
			int col = 3*face.Index();
			
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

			vcg::Point3f n01 = nabla_phi0^nabla_phi1;
			vcg::Point3f n12 = nabla_phi1^nabla_phi2;
			vcg::Point3f n20 = nabla_phi2^nabla_phi0;				

			float area = gdut_base::countArea(face);
			LL[std::make_pair(i0,i1)]+=n01*area;
			LL[std::make_pair(i1,i0)]-=n01*area;
			LL[std::make_pair(i1,i2)]+=n12*area;
			LL[std::make_pair(i2,i1)]-=n12*area;
			LL[std::make_pair(i2,i0)]+=n20*area;
			LL[std::make_pair(i0,i2)]-=n20*area;
		/*	LL.coeffRef(i0,i1)+=n01;
			LL.coeffRef(i1,i0)-=n01;

			LL.coeffRef(i0,i1)+=n01;
			L.coeffRef(i1,i0)-=n01;

			L.coeffRef(i0,i1)+=n01;
			L.coeffRef(i1,i0)-=n01;*/
		}
		L.setZero();
		for(unordered_map<std::pair<int, int>, vcg::Point3f>::iterator iter=LL.begin();iter!=LL.end();iter++){
			pair<int,int> indice = iter->first;
			vcg::Point3f p = iter->second;
			int row = indice.first*3;
			int col = indice.second*3;
			L.coeffRef(row,col+1) = -p.Z();
			L.coeffRef(row,col+2) = p.Y();
			L.coeffRef(row+1,col) = p.Z();
			L.coeffRef(row+1,col+2) = -p.Y();
			L.coeffRef(row+2,col) = -p.Y();
			L.coeffRef(row+2,col+1) = p.X();
		}
	
	//  L.makeCompressed();//ѹ��ʣ��Ŀռ�
}



void countDivfree(MeshModel &m,std::map<int,vcg::Point3f>& kexi,std::map<int,vcg::Point3f>& result){
	//    vector<float> x;// ��������|F|������
	int nV = m.cm.VN();
	int nF = m.cm.FN();
	VectorXd b(3*nV);
	b.setZero();
	buildRHS(m,kexi,b);

	double mean = b.mean();
	double maxb = b.maxCoeff();
	double minb = b.minCoeff();

	SparseMatrix<double> L(3*nV,3*nV);// ϵ������|3V|*|3V|��Ϊ��4�������Ĳ���
	L.reserve(36*nV);
	int nZ1 = L.nonZeros();
	buildLHS(m,L);
	int nZ = L.nonZeros();

	Eigen::SimplicialCholesky<SparseMatrix<double>> chol(L);  // ִ��A�� Cholesky�ֽ�
	//Eigen::ConjugateGradient<SparseMatrix<double>,Eigen::Upper> chol(A);  // ִ��A�� Cholesky�ֽ�
	VectorXd x(3*m.cm.VN());
	x = chol.solve(b);         // ʹ��A��Cholesky�ֽ������Ⱥ��ұߵ�����b

	double meanX = x.mean();
	double maxbX = x.maxCoeff();
	double minbX = x.minCoeff();

	if(chol.info()!=Success) {
		// decomposition failed
		for(int i = 0; i < m.cm.FN();i++){
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
			p0+=p1;
			p0+=p2;
			result.insert(pair<int,vcg::Point3f>(i,p0));
		}
		return;
	}	
}
}

#endif